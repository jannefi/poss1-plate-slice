#!/usr/bin/env bash
# Run on janne-pc ONLY, via cron (every 5 min) -- deliberately independent
# of any Claude session, SSH connection, or the laptop being open. Added
# after the 2026-08-25/26 incident: both full642 arms crashed when
# /dev/sdb2 filled to 0 bytes free (a capacity-planning error -- the tiles
# actually land there, not on /srv/vasco, which was the disk checked at
# planning time), and sat idle, billing, for ~12 hours before anyone
# noticed, because the only monitoring was a human manually asking for a
# status check. This script does not depend on that.
#
# Per arm: read its plate list + pulled-count to tell "finished" from
# "crashed"; if the orchestrator's own recorded PID is dead and its log has
# gone quiet, restart it (safe -- the whole pipeline is designed to resume
# idempotently via .pulled_<PLATE> markers and tile_status.json delta-skip).
# Backs off after repeated restarts with no progress rather than looping
# forever. Also checks disk space on both relevant mounts and each VM's
# gcloud-reported status.
#
# Install: crontab -l (check first) then add:
#   */5 * * * * /home/janne/code/poss1-plate-slice/scripts/gcp/watchdog.sh
# Cron runs with a minimal environment -- every external command below is
# called by absolute path, nothing relies on an inherited $PATH.
set -uo pipefail

GCLOUD=/usr/bin/gcloud
SSH=/usr/bin/ssh
DF=/usr/bin/df

REPO_ROOT="/home/janne/code/poss1-plate-slice"
KEY_PATH="/home/janne/.ssh/vasco60-gcp-pilot"
GCP_USER="janne"
PROJECT="project-a54d84d6-a5e7-4acc-b49"
ZONE="europe-north1-a"
VETO_WORKERS=3

LOG="$REPO_ROOT/work/runs/watchdog.log"
mkdir -p "$REPO_ROOT/work/runs"

NOW_TS="$(date -Is)"
NOTABLE=0
NOTES=()

note() {
  NOTABLE=1
  NOTES+=("$1")
  echo "[$NOW_TS] $1" >> "$LOG"
}

log() {
  echo "[$NOW_TS] $1" >> "$LOG"
}

# ---------------------------------------------------------------------
# Disk space, both relevant mounts.
# ---------------------------------------------------------------------
check_disk() {
  local PATH_TO_CHECK="$1" LABEL="$2" WARN_GB="$3" CRIT_GB="$4"
  local AVAIL_KB AVAIL_GB
  AVAIL_KB=$("$DF" --output=avail -k "$PATH_TO_CHECK" 2>/dev/null | tail -1 | tr -d ' ')
  if [ -z "$AVAIL_KB" ]; then
    note "disk check FAILED for $LABEL ($PATH_TO_CHECK) -- df returned nothing"
    return
  fi
  AVAIL_GB=$((AVAIL_KB / 1024 / 1024))
  if [ "$AVAIL_GB" -lt "$CRIT_GB" ]; then
    note "CRITICAL: $LABEL ($PATH_TO_CHECK) has only ${AVAIL_GB}GB free (< ${CRIT_GB}GB)"
  elif [ "$AVAIL_GB" -lt "$WARN_GB" ]; then
    note "WARNING: $LABEL ($PATH_TO_CHECK) has ${AVAIL_GB}GB free (< ${WARN_GB}GB)"
  else
    log "disk OK: $LABEL ${AVAIL_GB}GB free"
  fi
}

check_disk /srv/vasco "/srv/vasco (tiles target)" 500 100
check_disk / "/dev/sdb2 (OS/code)" 100 25

# ---------------------------------------------------------------------
# Per-arm liveness + auto-restart.
# ---------------------------------------------------------------------
check_arm() {
  local TAG="$1"
  local LOCAL_ROOT="$REPO_ROOT/work/runs/full642_gcp_${TAG}"
  local PLATES_FILE="$LOCAL_ROOT/plates.txt"
  local PIDFILE="$LOCAL_ROOT/orchestrator.pid"
  local HOSTFILE="$LOCAL_ROOT/gcp_host.txt"
  local INSTFILE="$LOCAL_ROOT/gcp_instance.txt"
  local ORCH_LOG="$LOCAL_ROOT/orchestrator.log"
  local RESTART_LOG="$LOCAL_ROOT/watchdog_restarts.log"

  [ -f "$PLATES_FILE" ] || { log "arm $TAG: no plates.txt, skipping (not started?)"; return; }

  local TOTAL DONE
  TOTAL=$(tr ',' '\n' < "$PLATES_FILE" | grep -c .)
  DONE=$(ls "$LOCAL_ROOT"/.pulled_* 2>/dev/null | wc -l)

  if [ "$DONE" -ge "$TOTAL" ]; then
    local SHUTDOWN_MARK="$LOCAL_ROOT/.vm_shutdown_done"
    if [ -f "$SHUTDOWN_MARK" ]; then
      log "arm $TAG: finished ($DONE/$TOTAL), VM already shut down -- nothing to do"
      return
    fi
    # VM's only job is steps 1-3 (slice); step4/5 runs entirely on janne-pc
    # using local catalog mirrors, so once every plate is pulled the VM has
    # nothing left to do, ever, for this arm -- safe to delete regardless of
    # step4/5 progress. Added 2026-08-27, explicitly approved by Janne:
    # previously this only logged "finished" and left the VM running (and
    # billing) until a human ran scripts/gcp/shutdown_vm.sh by hand -- risky
    # for an overnight finish while he's asleep.
    local INSTANCE=""
    [ -f "$INSTFILE" ] && INSTANCE=$(cat "$INSTFILE")
    if [ -z "$INSTANCE" ]; then
      note "arm $TAG: finished ($DONE/$TOTAL) but no $INSTFILE recorded -- cannot auto-shutdown, needs a human"
      return
    fi
    local STATUS
    STATUS=$("$GCLOUD" compute instances describe "$INSTANCE" --zone="$ZONE" --project="$PROJECT" \
      --format="value(status)" 2>/dev/null)
    if [ -z "$STATUS" ]; then
      log "arm $TAG: finished ($DONE/$TOTAL), VM $INSTANCE already gone"
      touch "$SHUTDOWN_MARK"
      return
    fi
    note "arm $TAG: finished ($DONE/$TOTAL) -- auto-deleting VM $INSTANCE to stop billing"
    if "$GCLOUD" compute instances delete "$INSTANCE" --zone="$ZONE" --project="$PROJECT" --quiet 2>>"$LOG"; then
      note "arm $TAG: VM $INSTANCE deleted"
      touch "$SHUTDOWN_MARK"
    else
      note "arm $TAG: VM $INSTANCE delete FAILED -- needs a human, see $LOG"
    fi
    return
  fi

  # VM state, informational -- doesn't gate the restart decision below,
  # since a stopped/preempted VM would just make the restart fail loudly
  # via SSH, which is itself useful signal in orchestrator.log.
  if [ -f "$INSTFILE" ]; then
    local INSTANCE STATUS
    INSTANCE=$(cat "$INSTFILE")
    STATUS=$("$GCLOUD" compute instances describe "$INSTANCE" --zone="$ZONE" --project="$PROJECT" \
      --format="value(status)" 2>/dev/null)
    if [ "$STATUS" != "RUNNING" ]; then
      note "arm $TAG: VM $INSTANCE status='$STATUS' (expected RUNNING)"
    fi
  fi

  local PID_ALIVE=0
  if [ -f "$PIDFILE" ]; then
    local PID; PID=$(cat "$PIDFILE")
    if [ -n "$PID" ] && kill -0 "$PID" 2>/dev/null; then
      PID_ALIVE=1
    fi
  fi

  if [ "$PID_ALIVE" -eq 1 ]; then
    # Alive is not the same as making progress. Added 2026-08-27 after a
    # router reboot left both orchestrators hung on a stale SSH channel --
    # the remote work had already finished (0% VM CPU) but the local ssh
    # client never got the completion signal, so it sat blocked forever
    # with its own PID still alive. The dead-pidfile check below never
    # fires for that case. Only a human noticing caught it that morning;
    # this closes the gap for when nobody's watching (e.g. overnight ISP
    # maintenance). Threshold checked live against a real slow-but-healthy
    # plate the same day: still going strong (fresh VM-side sex/psfex
    # processes, real progress) at 16+ minutes quiet log -- so this needs
    # real margin above "slow," not just above the ~5-15min typical case.
    # 2400s (40min) leaves that margin while still catching a genuine
    # all-night hang long before morning.
    local ORCH_LOG_AGE_S=999999
    if [ -f "$ORCH_LOG" ]; then
      ORCH_LOG_AGE_S=$(( $(date +%s) - $(stat -c %Y "$ORCH_LOG") ))
    fi
    if [ "$ORCH_LOG_AGE_S" -ge 2400 ]; then
      note "arm $TAG: STALLED (pid alive, $DONE/$TOTAL pulled, log quiet ${ORCH_LOG_AGE_S}s) -- killing so the next tick's dead-pidfile check can restart it cleanly"
      kill -9 "$PID" 2>/dev/null
      return
    fi
    log "arm $TAG: alive ($DONE/$TOTAL pulled)"
    return
  fi

  # Dead pidfile (or none). Avoid restarting something merely between
  # plates: only act if the log has been quiet for a while. ~13 min/plate
  # observed pace, with --progress writing continuously during a pull, so
  # 3 minutes of total silence should only happen when the process is
  # actually gone.
  local LOG_AGE_S=999999
  if [ -f "$ORCH_LOG" ]; then
    LOG_AGE_S=$(( $(date +%s) - $(stat -c %Y "$ORCH_LOG") ))
  fi
  if [ "$LOG_AGE_S" -lt 180 ]; then
    log "arm $TAG: pidfile dead but log active ${LOG_AGE_S}s ago -- not restarting yet"
    return
  fi

  # Backoff: 3+ restarts in the last 30 min with no net progress -> stop
  # auto-restarting, alert only. Prevents a silent crash-restart loop on
  # something genuinely broken from burning GCP compute unattended.
  local WINDOW_START=$(( $(date +%s) - 1800 ))
  local RECENT=0 FIRST_DONE=""
  if [ -f "$RESTART_LOG" ]; then
    while IFS=, read -r TS DONE_AT_RESTART; do
      if [ "${TS:-0}" -ge "$WINDOW_START" ] 2>/dev/null; then
        RECENT=$((RECENT+1))
        [ -z "$FIRST_DONE" ] && FIRST_DONE="$DONE_AT_RESTART"
      fi
    done < "$RESTART_LOG"
  fi
  if [ "$RECENT" -ge 3 ] && [ "$FIRST_DONE" = "$DONE" ]; then
    note "arm $TAG: BACKOFF -- $RECENT restarts in the last 30min with no progress (stuck at $DONE/$TOTAL). NOT auto-restarting again -- needs a human."
    return
  fi

  note "arm $TAG: CRASH DETECTED ($DONE/$TOTAL pulled, pidfile dead, log quiet ${LOG_AGE_S}s) -- auto-restarting"
  echo "$(date +%s),$DONE" >> "$RESTART_LOG"

  local GCP_HOST; GCP_HOST=$(cat "$HOSTFILE" 2>/dev/null)
  if [ -z "$GCP_HOST" ]; then
    note "arm $TAG: cannot restart -- no $HOSTFILE recorded"
    return
  fi
  local PLATES; PLATES=$(cat "$PLATES_FILE")
  local INSTANCE_ENV=""
  [ -f "$INSTFILE" ] && INSTANCE_ENV=$(cat "$INSTFILE")

  cd "$REPO_ROOT"
  nohup env GCP_HOST="$GCP_HOST" INSTANCE="$INSTANCE_ENV" PLATES="$PLATES" \
    LOCAL_ROOT="$LOCAL_ROOT" VETO_WORKERS="$VETO_WORKERS" KEY_PATH="$KEY_PATH" \
    /usr/bin/bash "$REPO_ROOT/scripts/gcp/run_and_pull_pilot.sh" \
    >> "$ORCH_LOG" 2>&1 < /dev/null &
  disown
  note "arm $TAG: restart launched (new pid $!)"
}

check_arm A
check_arm B

if [ "$NOTABLE" -eq 1 ]; then
  ALERT="$REPO_ROOT/work/runs/ALERT_$(date -u +%Y%m%dT%H%M%SZ).txt"
  {
    echo "watchdog alert at $NOW_TS"
    for n in "${NOTES[@]}"; do echo "- $n"; done
  } > "$ALERT"
fi

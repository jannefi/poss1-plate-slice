# Photographic Plate Forensics: a practical community library

> **Purpose.** Photographic survey plates are still scientifically valuable, but interpreting *digitised* plate products requires “plate literacy”: knowledge of PSF structure, photographic/emulsion artefacts, and digitisation/copying effects.
>
> This library collects **practical controls**, **diagnostic ideas**, and **hard-to-find references** so that discussions about plate-based anomalies can converge on reproducible checks rather than folklore.

## What this library is (and is not)
- **Is:** a curated set of references + reproducible “minimum viable controls” for evaluating plate-based anomalies.
- **Is not:** a claim about any single phenomenon; it aims to help *separate hypotheses* (optics/emulsion/digitisation vs astrophysical).

## Where the data come from (common sources)
- **Digitized Sky Survey (DSS):** digitised Schmidt survey plates distributed by STScI/MAST and mirrored elsewhere (DSS1/DSS2).  
- **APPLAUSE:** a modern digitisation + calibration archive with extracted source tables and metadata.  
- **MAPS/APS (POSS‑I):** digitised POSS‑I catalogues with computed image parameters and classification outputs.

## Minimum viable controls (MVC) for “anomalous” detections
1. **Cross-scanner / cross-digitisation provenance:** confirm presence/absence across independent scans (e.g., DSS vs SuperCOSMOS where available).
2. **Local PSF comparison:** compare candidate FWHM/shape/radial profile against *nearby stars of similar brightness*.
3. **Plate-position stratification:** check whether effects change with distance from plate centre / annular bins (focus/field effects).
4. **Edge buffers:** repeat analyses excluding plate edges/corners to reduce vignetting and processing artefact concentration.

## Common mechanism families (cheat sheet)
### A) Optics + PSF structure
- Gaussian-like cores plus non-Gaussian wings/scattered-light “aureoles”.
- Internal reflections, halation, ghost images.

### B) Photographic/emulsion + processing
- Development-related adjacency effects; local processing non-uniformities (“local errors”).
- Sensitivity to sky/background estimation at faint levels.

### C) Digitisation + copying procedures
- Scanner-dependent morphology differences.
- Copying/atlas reproduction artefacts (where relevant to the dataset).

## Diagnostics (quick plots and sanity checks)
- **Width vs shape departure:** candidate FWHM relative to local-star median vs a shape/departure metric.
- **Two-population check:** cluster in (width, shape) space, then inspect cutouts for each cluster.
- **Confidence regime awareness:** contour metrics on tiny areas can be unstable; stratify by confidence/area.

## Annotated bibliography (starter set)

### Core PSF / photographic-systematics references
- **Kormendy (1973)** — Star brightness profiles as standards; explicit discussion of ghost/reflection features and profile stability.  
  *Use for:* understanding Schmidt star-profile structure and instrumental artefacts. (DOI: https://ui.adsabs.harvard.edu/link_gateway/1973AJ.....78..255K/doi:10.1086/111412)

- **de Vaucouleurs (1984)** — Review of photographic photometry errors with Schmidt telescopes: PSF structure, internal reflections, scattered light (aureole), local errors, etc.  
  *Use for:* taxonomy of systematic effects and what can bias profile/shape measures. (Semantic Scholar:
https://www.semanticscholar.org/paper/Photometry-of-Extended-Sources-Vaucouleurs/4287cc9832072f089e454f91d64adb9efd5ef412)

- **King (1971)** - Core PSF / star-image formation (photographic plates & Schmidt systems).   
  *Use for:* Baseline star‑image radial profile on photographic material (core + exponential + inverse‑square aureole), including the reminder that the aureole carries ~5% of the light.
https://iopscience.iop.org/article/10.1086/129100

### Digitised survey context / access
- **DSS overview (MAST / ESO mirrors)** — what DSS is, what surveys it digitises, basic scanning characteristics.  
  *Use for:* dataset provenance and scan properties.
https://archive.stsci.edu/dss/

- **APPLAUSE documentation** — digitised plate archive with calibrations and extracted tables.  
  *Use for:* independent plate sources beyond POSS era/digitisation chain.
https://www.plate-archive.org/cms/home/

- **MAPS/APS (POSS‑I) overview** — digitised POSS‑I with computed image parameters and classification approaches.  
  *Use for:* morphology parameters at catalogue scale.
https://aps.umn.edu/docs/ 
https://iopscience.iop.org/article/10.1086/133186

### Calibration & Measurement Infrastructure
- **Humphreys, R. M., Landau, R., Ghigo, F. D., Zumach, W., & LaBonte, A. E. (1991)** BVR photoelectric sequences for selected fields in the Palomar Sky Survey and the magnitude–diameter relation. The Astronomical Journal, 102(1), 395+    
   *Use for* understanding how POSS‑I magnitudes and “image diameters” were calibrated using photoelectric BVR sequences and APS threshold‑densitometry. Documents plate‑to‑plate zero‑point shifts in the magnitude–diameter relation.
DOI: https://ui.adsabs.harvard.edu/link_gateway/1991AJ....102..395H/doi:10.1086/115883

### Defects in photographic plates
- **Greiner et al., Astronomy & Astrophysics, 234, 251–261 (1990)** Discrimination between star-like defects on photographic plates and possible gamma-ray burst optical counterparts.
    *Use for* Understanding plate defects.
DOI: https://ui.adsabs.harvard.edu/abs/1990A%26A...234..251G/abstract

### Modern morphology debates (examples of current work)
- **Hambly & Blair (2024)** — morphology-based classification using SuperCOSMOS; discusses possible copying-procedure artefacts.  
  *Use for:* cautionary reference on digitised “transient-like” sources.
https://academic.oup.com/rasti/article/3/1/73/7601398 

- **Busko (2026 preprint)** — independent archival-plate transient search; reports systematic narrow FWHM behaviour relative to stars in his dataset.  
  *Use for:* contemporary methodology and comparative diagnostics.
https://arxiv.org/abs/2603.20407 

## Contributing
- Please submit additions as **citations + short annotation**, not full PDFs.
- Prefer **DOIs / ADS bibcodes / arXiv**; avoid redistributing copyrighted full text unless explicitly permitted.
- If you add an “artefact mechanism”, include: dataset context, telltale signatures, and suggested controls.

## Licence / note
This repository contains *references and commentary*, not redistributed copyrighted publications.

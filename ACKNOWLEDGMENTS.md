# Acknowledgments and data provenance

> **These are conditions of use, not courtesies.** If you use these data or these
> derived products, the acknowledgments below are required.
>
> The wording follows each provider's standard text as of August 2026. Providers
> do revise it — check the current page before reusing these in a publication.

## Survey imagery

**Digitized Sky Survey (DSS).** Plate headers in this work record
`COPYRGHT = 'STScI/AURA'`.

> The Digitized Sky Surveys were produced at the Space Telescope Science
> Institute under U.S. Government grant NAG W-2166. The images of these surveys
> are based on photographic data obtained using the Oschin Schmidt Telescope on
> Palomar Mountain and the UK Schmidt Telescope. The plates were processed into
> the present compressed digital form with the permission of these institutions.

**POSS-I.**

> The National Geographic Society - Palomar Observatory Sky Atlas (POSS-I) was
> made by the California Institute of Technology with grants from the National
> Geographic Society.

**IRSA.** The full-plate scans are obtained from the NASA/IPAC Infrared Science
Archive, operated by the Jet Propulsion Laboratory, California Institute of
Technology, under contract with NASA.

**MAST / STScI.** Archive cutouts are obtained from the Mikulski Archive for
Space Telescopes.

## Reference catalogues

**Gaia.**

> This work has made use of data from the European Space Agency (ESA) mission
> Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing
> and Analysis Consortium (DPAC,
> https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has
> been provided by national institutions, in particular the institutions
> participating in the Gaia Multilateral Agreement.

**Pan-STARRS1.** The Pan-STARRS1 Surveys have been made possible through
contributions by the Institute for Astronomy, the University of Hawaii, and the
PS1 Science Consortium. Cite Chambers et al. (2016), arXiv:1612.05560.

**USNO-B1.0.** Monet, D. G., et al. 2003, AJ, 125, 984.

**VizieR / CDS.** This research has made use of the VizieR catalogue access tool
and the SIMBAD database, CDS, Strasbourg, France. Cite Ochsenbein et al. (2000),
A&AS, 143, 23.

## Software

| tool | citation |
|---|---|
| SExtractor | Bertin & Arnouts (1996), A&AS, 117, 393 |
| PSFEx | Bertin (2011), ASP Conf. Ser., 442, 435 |
| STILTS | Taylor (2006), ASP Conf. Ser., 351, 666 |
| Astropy | Astropy Collaboration (2013, 2018, 2022) |
| NumPy | Harris et al. (2020), Nature, 585, 357 |
| SciPy | Virtanen et al. (2020), Nature Methods, 17, 261 |
| pandas | McKinney (2010), Proc. 9th Python in Science Conf. |

## Astrometric model

The two-solution behaviour documented in `docs/DSS_WCS_TWO_SOLUTIONS.md` rests
on the GSSS plate model:

- Lasker, B. M., et al. 1990, AJ, 99, 2019
- Russell, J. L., et al. 1990, AJ, 99, 2059
- Morrison, J. E., et al. 2001, AJ, 121, 1752

## Relationship to prior work

This repository is an **independent** reproduction of a POSS-I vanishing-source
search. It is not affiliated with, endorsed by, or produced in collaboration with
the VASCO collaboration or any other group, and it uses **no** unpublished data
from any of them.

Validation is against the published POSS-I vanishing-source catalogue available
from SVO. Where results are compared with published work, that work is cited in
`RESULTS.md`.

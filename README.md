# Developmental and Phylogenetic Patterns of Antennal Expansion in Leaf-Footed Bugs (Hemiptera: Coreidae)

_Data and R code for the comparative study of antennal expansion in leaf-footed bugs, scored across adults and nymphs and reconstructed over a recent Coreidae phylogeny_

Authors: Mariana Polli, Ummat Somjee & Alexandre V. Palaoro <br>
Contact about code and analyses: marianapolli12@gmail.com or alexandre.palaoro@gmail.com

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22758792.svg)](https://doi.org/10.5281/zenodo.22758792)
[![preprint](https://img.shields.io/badge/preprint-10.32942%2FX21X1M-blue)](https://doi.org/10.32942/X21X1M)
[![code license](https://img.shields.io/badge/code%20license-MIT-green)](LICENSE)
[![data license](https://img.shields.io/badge/data%20license-CC0%201.0-brightgreen)](LICENSE-DATA)

> [!NOTE]
> This work is a preprint under peer review. This README will be updated with the journal
> reference once the paper is published.

---

### This readme has been divided in three parts. First, we will talk about file structure, then the code, then the dataset.

##### File structure:

We have four folders: <i>"code"</i>, <i>"data"</i>, <i>"evo.models"</i>, and <i>"phylo"</i>. <br>
The <i>"code"</i> folder contains one file, `coreid-analyses.R`, with all the code required to run the analyses. <br>
The <i>"data"</i> folder contains all the data required to run our analyses, in ".csv". All the data files begin with the last date they were updated, followed by a descriptive name. The file "Coreidae_Adults" contains all species data for adults. The files containing the word "genera" contain only one species per genus. The file containing the word "summarised" contains only the species we have information for and that were not excluded during our search. <br>
The <i>"evo.models"</i> folder contains saved evolutionary analyses. We saved them in different ".RData" files because some of them might take >5 minutes to run. Each file is a different analysis and they are called at different times in the code. <br>
The <i>"phylo"</i> folder contains the phylogeny we used in the paper. <br>

We also have four files that are not inside any folder:

- `Coreidae_Adults.csv` — contains the classification of each antennomere and links to the images we used to categorize antennal morphology in adults. Also has links to primary taxonomic descriptions in a few instances.
- `Coreidae_Nymphs_Full.csv` — contains the classification of each antennomere and links to the images we used to categorize antennomere morphology in nymphs. This is the full file; we collapsed it into genus for the analyses.
- `antennae-evol.Rmd` — the RMarkdown source of the report: the analyses as they were actually run, with comments that explains each step.
- `antennae-evol.html` — the knitted report produced from `antennae-evol.Rmd`.

##### Code:

There are two equivalent entry points: `antennae-evol.Rmd`, which knits the full report, and
`code/coreid-analyses.R`, the same analysis as a plain script. Both read the same data and load the
same stored model fits, and both use paths relative to the repository root — open
`coreidae-antennae` as the working directory (or as an RStudio project) and everything resolves.

The pipeline reads the phylogeny and the scored antennal data, roots the tree following Forthman et al. (2024), prunes it to the taxa with data, and fits equal-rates (ER), all-rates-different (ARD) and hidden-rates (HRM) Mk models under both equal and
FitzJohn root priors — separately for adults at species level, adults at genus level, and nymphs at
genus level. Model fits are stored in `evo.models/` so the script can be re-run without waiting for
every optimisation, and the final blocks build the stochastic-character maps and the paired
adult/nymph figure.

##### Dataset:

The datasheets containing Coreidae antennae classification share a similar format.

First column - Species ID in number; <br>
Second column - how the species is written at the tips of the phylogeny file, to match the information more easily; <br>
Third column - "Species_Name", mirroring how the name is written at the tips of the phylogeny file; <br>
Fourth column - Family, since some outgroups are not Coreidae; <br>
Fifth column - Subfamily; <br>
Sixth column - Tribe; <br>
Seventh column - Genus; <br>
Eighth column - Species; <br>
Ninth column - Antennae classification, whether we classified it as having expansion "2", no expansion "1", or no data "0"; <br>
Tenth column - State of antennomere I (varying from 0 to 2); <br>
Eleventh column - State of antennomere II (varying from 0 to 2); <br>
Twelfth column - State of antennomere III (varying from 0 to 2); <br>
Thirteenth column - Sum of antennomeres (varying 0 to 6); <br>
Fourteenth column - Reference of images or papers used to classify the species' antennae. <br>

The main difference is that in the full datasets, antennae can also be classified as "No_data". <br>

The CSV files use **semicolon** as the field separator (`sep = ";"`), not comma.

##### Phylogeny:

> Forthman, M., Phan, H., Miller, C. W. & Kimball, R. T. (2024) Phylogenetic placement of the leaf-footed bug tribes Agriopocorini, Amorbini, and Manocoreini (Heteroptera: Coreidae) using ultraconserved elements. *Zoological Journal of the Linnean Society* 202(3): zlae024. [https://doi.org/10.1093/zoolinnean/zlae024](https://doi.org/10.1093/zoolinnean/zlae024)

##### Packages:

The code was run in R software v4.6.0. <br>
Packages used: <br>
- pander (v0.6.6) <br>
- lubridate (v1.9.5) <br>
- forcats (v1.0.1) <br>
- stringr (v1.6.0) <br>
- dplyr (v1.2.1) <br>
- purrr (v1.2.2) <br>
- readr (v2.2.0) <br>
- tidyr (v1.3.2) <br>
- tibble (v3.3.1) <br>
- ggplot2 (v4.0.3) <br>
- tidyverse (v2.0.0) <br>
- diversitree (v0.10-1) <br>
- geiger (v2.0.11) <br>
- phytools (v2.5-2) <br>
- maps (v3.4.3) <br>
- ape (v5.8-1)

---

## Citation

If you use anything in this repository, please cite the paper:

> Polli, M., Somjee, U. & Palaoro, A.V. (2026) Developmental and phylogenetic patterns of antennal expansion in leaf-footed bugs (Hemiptera: Coreidae). *EcoEvoRxiv*. [https://doi.org/10.32942/X21X1M](https://doi.org/10.32942/X21X1M)

This is a preprint under peer review — please check for a journal version before citing.

If you reuse the code or the archived files directly, please also cite the archive:

> Polli, M., Somjee, U. & Palaoro, A.V. (2026) *Code and data for: Developmental and Phylogenetic Patterns of Antennal Expansion in Leaf-Footed Bugs (Hemiptera: Coreidae)* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.22758792

`CITATION.cff` in this repository holds the citation in machine-readable form — GitHub's
**Cite this repository** button (top right of the repository page) will generate APA or
BibTeX from it for you.

## License

This repository is released under two licenses, because code and data are different things.

| Content | License | File |
| --- | --- | --- |
| Analysis code — `code/coreid-analyses.R` and the knitted `antennae-evol.html` | [MIT](https://opensource.org/licenses/MIT) | [`LICENSE`](LICENSE) |
| Data, phylogeny and stored model output — everything in `data/`, `phylo/` and `evo.models/`, plus the two root-level `.csv` files | [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) | [`LICENSE-DATA`](LICENSE-DATA) |

In short: do whatever you like with the data, no permission needed and no attribution
legally required; reuse the code freely as long as you keep the copyright notice. Academic
norms still apply — if the data or code are useful to you, cite the paper.

## Reproducibility

- Everything needed to reproduce the analyses is in this repository: the scored data, the
  phylogeny, the stored model fits and the script that produced them.
- The archived release on Zenodo is the version of record. GitHub history may move on; the
  DOI will not.
- Package versions used are listed above. If a result does not reproduce, check those first.
- The stochastic character maps (`evo.models/ant_sim.RData`, `evo.models/ard_sim_nymph.RData`) are
  the 1000 simulations behind the published figures. The `make.simmap()` calls that generated them
  are kept in the code but commented out, since re-running them is slow and gives a different
  draw each time; uncomment them if you want a fresh set.
- Likewise, the `save()` calls that wrote the files in `evo.models/` are commented, so running
  the code cannot overwrite the archived fits.
- Figures are plotted to the screen. To write them to file instead, uncomment the `pdf()`/`png()`
  and matching `dev.off()` lines and create a `figures/` folder first.

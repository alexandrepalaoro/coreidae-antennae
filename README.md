# Evolutionary history suggests nymphal origin of antennal expansions in adult leaf-footed bugs

_Data and R code for the comparative study of antennal expansion in leaf-footed bugs (Hemiptera: Coreidae), scored in adults and nymphs and reconstructed over a recent Coreidae phylogeny._

Authors: Mariana Polli, Ummat Somjee & Alexandre V. Palaoro <br>
Contact about code and analyses: marianapolli12@gmail.com or alexandre.palaoro@gmail.com

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22758792.svg)](https://doi.org/10.5281/zenodo.22758792)
[![article](https://img.shields.io/badge/article-10.1111%2Feen.70159-blue)](https://doi.org/10.1111/een.70159)
[![code license](https://img.shields.io/badge/code%20license-MIT-green)](LICENSE)
[![data license](https://img.shields.io/badge/data%20license-CC0%201.0-brightgreen)](LICENSE-DATA)

## The paper

> Polli, M., Somjee, U. & Palaoro, A.V. (2026) Evolutionary history suggests nymphal origin of antennal expansions in adult leaf-footed bugs. *Ecological Entomology*. [https://doi.org/10.1111/een.70159](https://doi.org/10.1111/een.70159)

The article is open access. In short: expanded antennae are rare in adult Coreidae and restricted to Coreinae in our sample, yet they arose several times. Nymphs show expansion far more often, including in genera whose adults have plain, rectilinear antennae. Every genus with expanded adult antennae also has expanded nymphal antennae.

This repository holds everything needed to repeat the analyses behind that result.

---

### This readme has three parts: file structure, then the code, then the dataset.

##### File structure:

We have four folders: <i>"code"</i>, <i>"data"</i>, <i>"evo.models"</i>, and <i>"phylo"</i>. <br>
The <i>"code"</i> folder contains one file, `coreid-analyses.R`, with all the code required to run the analyses. <br>
The <i>"data"</i> folder contains the data required to run our analyses, in ".csv". Each file name begins with the date of its last update, followed by a descriptive name. `Coreidae_Adults` holds all adult species. The files with "genera" in the name hold one species per genus. The file with "summarised" in the name lists the species on the phylogeny, in tip order, so the tip labels can be swapped for species names. <br>
The <i>"evo.models"</i> folder contains saved evolutionary analyses. We saved them in separate ".RData" files because some take more than 5 minutes to run. Each file is a different analysis, and the code loads them where needed. <br>
The <i>"phylo"</i> folder contains the phylogeny we used in the paper. <br>

We also have four files that are not inside any folder:

- `Coreidae_Adults.csv`: the classification of each antennomere in adults, with links to the images we used. In a few cases it also links to primary taxonomic descriptions.
- `Coreidae_Nymphs_Full.csv`: the same for nymphs. This is the full file; we collapsed it to genus for the analyses.
- `antennae-evol.Rmd`: the RMarkdown source of the report. It has the analyses as we ran them, with the prose that explains each step.
- `antennae-evol.html`: the knitted report produced from `antennae-evol.Rmd`.

There is also `antennae-evol.Rproj`, the RStudio project file. Open it and every path below resolves.

##### Code:

There are two equivalent entry points: `antennae-evol.Rmd`, which knits the full report, and
`code/coreid-analyses.R`, the same analysis as a plain script. Both read the same data and load the
same stored model fits. Both use paths relative to the repository root, so open `coreidae-antennae`
as the working directory (or as an RStudio project) before you run anything.

The pipeline reads the phylogeny and the scored antennal data, roots the tree following Forthman et al. (2024), prunes it to the taxa with data, and fits
equal-rates (ER), all-rates-different (ARD) and hidden-rates (HRM) Mk models under both an equal and a
FitzJohn root prior. It does this separately for adults at species level, adults at genus level, and nymphs at
genus level. Models are ranked by AIC and Akaike weight. The final blocks build the stochastic character maps and the paired
adult/nymph figure.

The script returns the AIC tables reported in the paper (Tables S1 to S3):

| Dataset | Best model | AIC | Akaike weight | Paper |
| --- | --- | --- | --- | --- |
| Adults, species level (174 species) | ER, FitzJohn root | 88.16 | 0.61 | Table S1 |
| Nymphs, genus level (87 genera) | HRM-ER, FitzJohn root | 137.95 | 0.89 | Table S2 |
| Adults, genus level | ER, FitzJohn root | 71.96 | 0.61 | Table S3 |

##### Dataset:

The datasheets with the antennae classification share a similar format. Columns, in order:

- First column: a row ID joined to the tip label of the phylogeny. It loads as a column named `X..x.`;
- `Species_Name`: how the species is written at the tips of the phylogeny file, so the data match the tree;
- `Family`: Coreidae, or the family of an outgroup;
- `Subfamily`;
- `Tribe`;
- `Genus`;
- `Species`;
- `Instar` (nymph files only): the lowest instar recorded for the genus. Antennal expansion can change across nymphal stages, so we used the earliest one;
- `Antennae`: the overall classification. `1` is no expansion, `2` is at least one antennomere expanded, `0` is no data;
- `Antennomere_I`, `Antennomere_II`, `Antennomere_III`: the state of each antennomere. `1` is rectilinear, `2` is expanded, `0` is no data. We did not score antennomere IV because it did not vary;
- `Sum`: the sum of the three antennomeres. It runs from 3 (none expanded) to 6 (all three expanded), and it is the character state used in the models. Expansion was nested: antennomere III is expanded whenever any is, II is added when two are, and I only when all three are, so the sum identifies the combination;
- `Reference`: the iNaturalist observation or the paper used to classify the antennae.

Species with no usable image or description, and tips known only to genus, have `0` and were dropped before the analyses.

The CSV files in `data/` use **semicolons** as the field separator for `18-08-26_Coreidae_Adults.csv` and `25-10-25_antena_coreidae_summarised_2.csv`, and **commas** for the two `*_genera.csv` files. The code reads each one with the right separator.

##### Phylogeny:

> Forthman, M., Phan, H., Miller, C. W. & Kimball, R. T. (2024) Phylogenetic placement of the leaf-footed bug tribes Agriopocorini, Amorbini, and Manocoreini (Heteroptera: Coreidae) using ultraconserved elements. *Zoological Journal of the Linnean Society* 202(3): zlae024. [https://doi.org/10.1093/zoolinnean/zlae024](https://doi.org/10.1093/zoolinnean/zlae024)

##### Packages:

The report was knitted in R v4.6.1. <br>
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
- geiger (v2.0.12) <br>
- phytools (v2.5-2) <br>
- maps (v3.4.3) <br>
- ape (v5.8-1)

---

## Citation

If you use anything in this repository, please cite the paper:

> Polli, M., Somjee, U. & Palaoro, A.V. (2026) Evolutionary history suggests nymphal origin of antennal expansions in adult leaf-footed bugs. *Ecological Entomology*. [https://doi.org/10.1111/een.70159](https://doi.org/10.1111/een.70159)

If you reuse the code or the archived files directly, please also cite the archive:

> Polli, M., Somjee, U. & Palaoro, A.V. (2026) *Code and data for: Evolutionary history suggests nymphal origin of antennal expansions in adult leaf-footed bugs* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.22758792

An earlier version of the paper was posted as a preprint ([EcoEvoRxiv, 10.32942/X21X1M](https://doi.org/10.32942/X21X1M)). Cite the journal article instead.

`CITATION.cff` in this repository holds the citation in machine-readable form. GitHub's
**Cite this repository** button (top right of the repository page) will generate APA or
BibTeX from it for you.

## License

This repository is released under two licenses, because code and data are different things.

| Content | License | File |
| --- | --- | --- |
| Analysis code: `code/coreid-analyses.R`, `antennae-evol.Rmd` and the knitted `antennae-evol.html` | [MIT](https://opensource.org/licenses/MIT) | [`LICENSE`](LICENSE) |
| Data, phylogeny and stored model output: everything in `data/`, `phylo/` and `evo.models/`, plus the two root-level `.csv` files | [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) | [`LICENSE-DATA`](LICENSE-DATA) |

In short: do whatever you like with the data, no permission needed and no attribution
legally required. Reuse the code freely as long as you keep the copyright notice. Academic
norms still apply: if the data or code are useful to you, cite the paper.

## Updating the archive

This repository is linked to Zenodo. Every new GitHub release is archived automatically and
gets its own version DOI. The concept DOI above always resolves to the newest version, so the
DOI in this README and in `CITATION.cff` never needs changing.

To publish an update: **Releases → Draft a new release**, bump the tag, publish.
Zenodo picks it up within a minute or two. `.zenodo.json` supplies the title, authors, ORCIDs,
keywords and the links to the article and preprint, so there is nothing to retype in the Zenodo form.

Cite the concept DOI in papers, never a version DOI.

## Reproducibility

- Everything needed to reproduce the analyses is in this repository: the scored data, the
  phylogeny, the stored model fits and the script that produced them.
- The script runs from start to finish on a clean copy of the repository and returns the AIC tables in the paper (table above).
- The archived release on Zenodo is the version of record. GitHub history may move on; the
  DOI will not.
- Package versions are listed above. If a result does not reproduce, check those first.
- The stochastic character maps (`evo.models/ant_sim.RData` for adults, `evo.models/hrm_sim_nymph.RData` for nymphs) are
  the 1000 simulations behind the published figures. The `make.simmap()` and `simmap()` calls that generated them
  are kept in the code but commented out, since re-running them is slow and gives a different
  draw each time. Uncomment them if you want a fresh set.
- Likewise, the `save()` calls that wrote the files in `evo.models/` are commented out, so running
  the code cannot overwrite the archived fits.
- Figures are plotted to the screen. To write them to file instead, uncomment the `pdf()`/`png()`
  and matching `dev.off()` lines and create a `figures/` folder first.

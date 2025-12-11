# Evolution and Phylogenetic Conservatism of Antennal Morphology in Leaf-Footed Bugs (Coreidae)

Authors: Mariana Polli, Ummat Somjee & Alexandre V. Palaoro <br>
---

Contact about code and analyses: marianapolli@ufpr.br or alexandre.palaoro@gmail.com <br>

---

This paper is under evaluation in a scientific journal. If you use any of our code or data, please email us just so we are aware <br>


### This readme has been divided in three parts. First, we will talk about file structure, then the code, the dataset.

##### File structure:

We have three folders: <i>"code"</i>, <i>"data"</i>, <i>"evo.models"</i>, and "phylo".
The <i>"code"</i> folder contains one file with all the code required to run the analyses. <br>
The <i>"data"</i> folder contains all the data required to run our analyses. In it, you will find the datasets in ".csv". All the data files begin with the last date they were updated, and a descriptive name for the file. The file "Coreidae antena" contains all species in the phylogeny. The file containing the word "summarized" contain only the species we have information for and were not excluded during our search. The file containing the word "uncertain" contains the species we also do not have information for. <br>
The <i>"evo.models"</i> folder contains saved evolutionary analyses. We saved them in different ".Rdata" files because some of them might take >5 minutes to run. Each file is a different analysis and they are called at different times in the code. <br>
The <i>"phylo"</i> folder contains the phylogeny we used in the paper.
We also have four files that are not inside any folder.
- "All_records_Coreidae_antena.csv" - Contains all the species with the links to the images we used to categorize the species.
- "antennae-evol.html" - an Rmarkdown file with all the analyses performed.
- "corr-t2-antena.Rdata" - the fossil-calibrated phylogeny using the "correlated" method from the package ape.
- "relax-t2-antena.Rdata" - the fossil-calibrated phylogeny using the "relaxed" method from the package ape.


##### Dataset:

The datasheets containing Coreidae antennae classification contain similar formats.

- First column - Species ID in number; <br>
- Second column - how the species is written is the tips of the phylogeny file to match the information more easily;
- Third column - "Species_name" that mirror how the name is written in the tips of the phylogeny file;
- Fourth column - Family, since some outgroups are not Coreidae; <br>
- Fifth column - Subfamily; <br>
- Sixth column - Tribe; <br>
- Seventh column - Genus; <br>
- Eighth column - Species; <br>
- Ninth column - Antennae classification, whether we classified it as having expansion "Expansion", or having no expansion "Straight" <br>

The main difference is that in the uncertain and full datasets, antennae can also be classified as "No_data". <br>

For the calibration points, we used the format suggested by the package "ape", where:
- First column - "node", the number of the node in which the age should be place (i.e., where the fossil should be placed);
- Second column - "age.min", the minimum age for that node;
- Third column - "age.max", the maximum age for that node;
- Fourth column - "soft.bounds", unused but required to run according to the function vignette.
For information on which fossils we used, please check the Supplementary Information. 

## Reference

Forthman, M., Phan, H., Miller, C. W., & Kimball, R. T. (2024). Phylogenetic placement of the leaf-footed bug tribes Agriopocorini, Amorbini, and Manocoreini (Heteroptera: Coreidae) using ultraconserved elements. Zoological Journal of the Linnean Society, 202(3), zlae024. 

## Packages 

The code was run in R software v4.5.1. <br>
Packages used: <br>
- lubridate(v.1.9.4) <br>
- forcats(v.1.0.1) <br>
- stringr(v.1.5.2) <br>
- dplyr(v.1.1.4) <br>
- purrr(v.1.1.0) <br>
- readr(v.2.1.5) <br>
- tidyr(v.1.3.1) <br>
- tibble(v.3.3.0) <br>
- ggplot2(v.4.0.0) <br>
- tidyverse(v.2.0.0) <br>
- pander(v.0.6.6) <br>
- diversitree(v.0.10-1) <br>
- geiger(v.2.0.11) <br> 
- phytools(v.2.4-4) <br>
- maps(v.3.4.2.1) <br>
- ape(v.5.8-1)

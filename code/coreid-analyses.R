library(ape)
library(phytools)
library(geiger)
library(diversitree)
library(tidyverse)

#### LOADING PHYLOGENY AND DATA FILE ####

tt = read.tree(file = './phylo/50p.treefile')

tips = read.csv(file = './data/antena_coreidae_summarised_2.csv', h =T, sep = ';')

#### LOADING  ADULTS (SPECIES), NYMPHS (GENUS) AND ADULTS (GENUS) DATA ####
sp_adults = read.csv("./data/Coreidae_Adults.csv", sep = ";") %>% 
  mutate(Species_Name = case_when(
    Species_Name == "Thasus_neocalifornicus" ~ "Pachylis_neocalifornicus",
    Species_Name == "Thasus_sp." ~ 'Pachylis_sp.',
    TRUE ~ Species_Name
  ))

gen_adults = read.csv("./data/Adults_genera.csv", sep = ";") %>% 
  mutate(Species_Name = case_when(
    Species_Name == "Thasus_neocalifornicus" ~ "Pachylis_neocalifornicus",
    Species_Name == "Thasus_sp." ~ 'Pachylis_sp.',
    TRUE ~ Species_Name
  ))

gen_nymphs = read.csv('./data/Nymphs_genera.csv', sep = ';') %>% 
  mutate(Species_Name = case_when(
    Species_Name == "Thasus_neocalifornicus" ~ "Pachylis_neocalifornicus",
    Species_Name == "Thasus_sp." ~ 'Pachylis_sp.',
    TRUE ~ Species_Name
  ))

#### HANDLING THE NAMES ON THE TIPS OF THE TREE AND ON OUR DATA FILE ####

rownames(tips) = tips$Species_Name

tt$tip.label = tips$Species_Name

name.check(tt, tips)

#### ROOTING THE TREE ####

#### TO GET A RESULT CLOSER TO WHAT WAS PRESENTED IN THE PAPER, 
#### AND TO BE ABLE TO DATE THE TREE, IT FIRST NEEDS A ROOT 
#### FOLLOWING THE INFORMATION ON FORTHMANN ET AL. (2024), 
#### JADERA IS THE ROOT. SO, LET'S ROOT IT.

t1 = root(tt, node = 308, resolve.root = T)

plot.phylo(ladderize(t1))
plot.phylo(t1)

#### ADULTS ANALYSES AT SPECIES LEVEL ####

#### REMOVING ZEROES FROM DATASET AND DROPPING UNNECESSARY TIPS ####

# Since our goal is analyzing trait evolution in Coreidae,
# we removed the non-Coreidae species in the tree.

# We also removed tips with no data or at species level,
# labeled as "0"

adults_data = sp_adults %>% 
  filter(Antennae != 0, 
         !Species_Name %in% c("Jadera_haematoloma", "Harmostes_serratus",
                              "Oedancala_sp.", "Oncopeltus_fasciatus",
                              "Largus_sp.", "Dysdercus_mimus",
                              "Dysdercus_suturellus", "Halyomorpha_halys"))


# Now we prune the tree to the data we have on Coreidae adults.

t_adults = keep.tip(t1, adults_data$Species_Name)

name.check(t_adults, adults_data)

#### THE MK MODEL: ASSESSING CHANGES IN ANTENNOMMERE ####
#### TRAIT STATE ON ADULTS SPECIES ####

# Before running the analyses, we need to create a matrix so phytools
# can read our data.


adults.sp.sim = setNames(as.factor(adults_data$Sum), adults_data$Species_Name)

adults.sp.sim

# With the matrix done, we can now move to the analyses.

# We are using the fitMk function. The function requires us to specify 
# the tree we are going to use, the data in matrix format, the evolutionary
# model, and, if we wish, the state of the root.
# We opted to calculate the transitions for the tree calibrated with the 
# correlated model in chronos() because it is most comonnly used model.

# We are fitting two types of evolutionary models: Equal rates and All rates differ.
# We are not using the symmetrical because we only have two trait states. 
# When that happens, the symmetrical evolutionary model is equal to the all rates differ
# model.

# Lastly, we are using two types of root estimation. One where both states have
# equal probabilities, and another where probabilities are calculated based on the
# data (fitzjohn). We used both because we have no a priori hypothesis of the
# ancestral state.

### NOTE: ONCE AGAIN I AM COMMENTING THE FUNCTION SO YOU DO NOT MISTAKENLY
### RUN IT. IF YOU WANT TO RUN IT, JUST UNCOMMENT THE FUNCTIONS. 
### OTHERWISE, LOAD THEM AND CARRY ON WITH THE ANALYSES.

#fit.ER.equal = fitMk(ladderize(t_adults), adults.sp.sim , model = "ER",
#                     pi = 'equal')


#fit.ER.fitz = fitMk(ladderize(t_adults), adults.sp.sim , model = "ER",
#                    pi = 'fitzjohn')

#plot(fit.ER.fitz)

#save(fit.ER.equal, file = './evo.models/er-equal-antenomero.RData')
#save(fit.ER.fitz, file = './evo.models/er-fitz-antenomero.RData')

load( file = './evo.models/er-equal-antenomero.RData')
load( file = './evo.models/er-fitz-antenomero.RData')

plot(fit.ER.equal)
fit.ER.equal

plot(fit.ER.fitz)
fit.ER.fitz

# All rates differ model #

#fit.ARD.equal = fitMk(ladderize(t_adults), adults.sp.sim., model = "ARD",
#                      pi = 'equal')

#fit.ARD.fitz = fitMk(ladderize(t_adults), adults.sp.sim, model = "ARD",
#                     pi = 'fitzjohn')

#save(fit.ARD.equal, file = './evo.models/ard-equal-antenomero.RData')
#save(fit.ARD.fitz, file = './evo.models/ard-fitz-antenomero.RData')

load(file = './evo.models/ard-equal-antenomero.RData')
load(file = './evo.models/ard-fitz-antenomero.RData')

plot(fit.ARD.equal)
fit.ARD.equal

plot(fit.ARD.fitz)
fit.ARD.fitz

# We are also using the HRM models for thoroughness. Maybe there is an
# underlying mechanism for changing antennal state we are not capturing. 

#fitHRM.er.equal = fitHRM(ladderize(t_adults), adults.sp.sim , model = "ER",
#                         pi = 'equal')

#fitHRM.er.fitz = fitHRM(ladderize(t_adults), adults.sp.sim , model = "ER",
#                        pi = 'fitzjohn')

#save(fitHRM.er.equal, file = './evo.models/hrm-er-equal-antenomero.RData')
#save(fitHRM.er.fitz, file = './evo.models/hrm-er-fitz-antenomero.RData')

load(file = './evo.models/hrm-er-equal-antenomero.RData')
load(file = './evo.models/hrm-er-fitz-antenomero.RData')

plot(fitHRM.er.equal)
fitHRM.er.equal

plot(fitHRM.er.fitz)
fitHRM.er.fitz

# All rates differ HRM #

#fitHRM.equal = fitHRM(ladderize(t_adults), adults.sp.sim , model = "ARD",
#                      pi = 'equal')


#fitHRM.fitz = fitHRM(ladderize(t_adults), adults.sp.sim , model = "ARD",
#                     pi = 'fitzjohn')

#save(fitHRM.equal, file = './evo.models/hrm-equal-antenomero.RData')
#save(fitHRM.fitz, file = './evo.models/hrm-fitz-antenomero.RData')


load(file = './evo.models/hrm-equal-antenomero.RData')
load(file = './evo.models/hrm-fitz-antenomero.RData')


plot(fitHRM.equal)
fitHRM.equal

plot(fitHRM.fitz)
fitHRM.fitz

# Now we compare the models using AIC. 
# We chose the model with the lowest AIC as the best model.
# If models had similar AICs, we chose the one with the lowest
# degrees of freedom. In the case, below, the ER model was the best.
# Between the two root states, we used the fitjohn approach because it 
# had more AIC weight when compared to the equal ER model.

fit_models = anova(fit.ER.equal, fit.ARD.equal,
                    fitHRM.er.equal, fitHRM.equal,
                    fit.ER.fitz, fit.ARD.fitz,
                    fitHRM.er.fitz, fitHRM.fitz)
  
anc.ARD = ancr(fit.ER.fitz)

fit.ER.fitz

cols = setNames(c("#0E3768","#94129F", "#de627a", "#ce4326"),
               levels(adults.sp.sim))


pdf(file = './figures/ancestral-tree-adults-species.pdf', h = 24, w = 20)
plot(anc.ant,  
     args.plotTree = list(lwd = 3, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          cex = 0.2, offset = 1
     ),
     args.nodelabels = list(cex = 0.3, piecol = cols),
     args.tiplabels = list(cex = 0.3)
     
)
dev.off()

# simulation mapping

ard_sim = make.simmap( tree = ladderize(t_adults), x= adults.sp.sim , model = "ER", pi = 'fitzjohn',
                     nsim = 1000)

sum_sim = describe.simmap(ard_sim)
sum_sim

plot(density(ard_sim))


#### NYMPHS AND ADULTS GENNERA ANALYSIS ####

# The same process made to antennomere analysis made for adults species will
# be made for nymphs and adults at the genus level

#### NYMPHS ####

# filtering data, removing 0 and species out of Coreidae

nymph_data = gen_nymphs%>% 
  filter(Antennae != 0, 
         !Species_Name %in% c("Jadera_haematoloma", "Harmostes_serratus",
                              "Oedancala_sp.", "Oncopeltus_fasciatus",
                              "Largus_sp.", "Dysdercus_mimus",
                              "Dysdercus_suturellus", "Halyomorpha_halys")) %>% 
  mutate(Species_Name = case_when(
    Species_Name == "Holhymenia_tibialis" ~ "Holhymenia_tibialis / clavigera",
    TRUE ~ Species_Name
  ))

# Now we prune the tree to the data we have on Coreidae nymphs

t_nymphs = keep.tip(t1, nymph_data$Species_Name)

name.check(t_nymphs, nymph_data)


# Creating nymphs matrix so phytools can read our data

nymphs.gen.sim = setNames(as.factor(nymph_data$Sum), nymph_data$Species_Name)

nymphs.gen.sim

# Model Mk, all models

#fit.ER.equal.nymph = fitMk(ladderize(t_nymphs), nymphs.gen.sim, model = "ER",
#                   pi = 'equal')


#fit.ER.fitz.nymph= fitMk(ladderize(t_nymphs), nymphs.gen.sim , model = "ER",
#                    pi = 'fitzjohn')

fit.ER.equal.nymph
fit.ER.fitz.nymph

plot(fit.ER.equal.nymph)

save(fit.ER.equal.nymph, file = './evo.models/er-equal-nymph.RData')
save(fit.ER.fitz.nymph, file = './evo.models/er-fitz-nymph.RData')

load( file = './evo.models/er-equal-nymph.RData')
load( file = './evo.models/er-fitz-nymph.RData')

#fit.ARD.equal.nymph = fitMk(ladderize(t_nymphs), nymphs.gen.sim , model = "ARD",
#                            pi = 'equal')

#fit.ARD.fitz.nymph = fitMk(ladderize(t_nymphs), nymphs.gen.sim, model = "ARD",
#                           pi = 'fitzjohn')

save(fit.ARD.equal.nymph, file = './evo.models/ARD-equal-nymph.RData')
save(fit.ARD.fitz.nymph, file = './evo.models//ARD-fitz-nymph.RData')

load( file = './evo.models/ARD-equal-nymph.RData')
load( file = './evo.models/ARD-fitz-nymph.RData')

# Hidden Rates Model, HRM

#fitHRM.er.equal.nymph = fitHRM(ladderize(t_nymphs), nymphs.gen.sim , model = "ER",
#                       pi = 'equal')

#fitHRM.er.fitz.nymph = fitHRM(ladderize(t_nymphs), nymphs.gen.sim , model = "ER",
#                              pi = 'fitzjohn')

save(fitHRM.er.equal.nymph, file = './evo.models/hrm-er-equal-nymph.RData')
save(fitHRM.er.fitz.nymph, file = './evo.models/hrm-er-fitz-nymph.RData')

load(file = './evo.models/hrm-er-equal-nymph.RData')
load(file = './evo.models/hrm-er-fitz-nymph.RData')

#fitHRM.ard.equal.nymph = fitHRM(ladderize(t_nymphs), nymphs.gen.sim , model = "ARD",
#                                pi = 'equal')


#fitHRM.ard.fitz.nymph = fitHRM(ladderize(t_nymphs), nymphs.gen.sim , model = "ARD",
#                         pi = 'fitzjohn')

save(fitHRM.ard.equal.nymph, file = './evo.models/hrm-equal-nymph.RData')
save(fitHRM.ard.fitz.nymph, file = './evo.models/hrm-fitz-nymph.RData')


load(file = './evo.models/hrm-equal-nymph.RData')
load(file = './evo.models/hrm-fitz-nymph.RData')

# comparing models
fit_models_nymphs = anova(fit.ER.equal.nymph, fit.ARD.equal.nymph,
                       fitHRM.er.equal.nymph, fitHRM.ard.equal.nymph,
                       fit.ER.fitz.nymph, fit.ARD.fitz.nymph,
                       fitHRM.er.fitz.nymph, fitHRM.ard.fitz.nymph)

# The chosen method was fit.HRM.er.fitz.nymph once it had the lowe AIC and
# less degrees of freedom (d.f. = 1)

fitHRM.er.fitz.nymph
plot(fitHRM.er.fitz.nymph)

anc.nymph = ancr(fitHRM.er.fitz.nymph)

cols = setNames(c("#0E3768", "#DE627A", "green", "brown"),
                levels(nymphs.gen.sim))


pdf(file = './ancestral-nymph.pdf', h = 40, w = 20)
plot(anc.nymph,  
     args.plotTree = list(lwd = 3, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          cex = 0.2, offset = 3
     ),
     args.nodelabels = list(cex = 0.3, piecol = cols),
     args.tiplabels = list(cex = 0.3)
     
)
dev.off()

# simulation mapping

ard_sim_nymph = make.simmap( tree = ladderize(t_nymphs), x= nymphs.gen.sim , model = "ER", pi = 'fitzjohn',
                             nsim = 1000)

sum_sim_nymph = describe.simmap(ard_sim_nymph)
sum_sim_nymph

dev.off()

plot(density(ard_sim_nymph))


#### ADULTS - GENUS LEVEL ####

# Now, we will use the adults data colapsed at genus level

adults_data2 = gen_adults %>% 
  filter(Antennae != 0, 
         !Species_Name %in% c("Jadera_haematoloma", "Harmostes_serratus",
                              "Oedancala_sp.", "Oncopeltus_fasciatus",
                              "Largus_sp.", "Dysdercus_mimus",
                              "Dysdercus_suturellus", "Halyomorpha_halys")) %>% 
  mutate(Species_Name = case_when(
    Species_Name == "Holhymenia_tibialis" ~ "Holhymenia_tibialis / clavigera",
    TRUE ~ Species_Name
  ))

# We need to make sure that only species present in nymphal data are present
# in the adults data, so we can compare them equally

adults_data2 = gen_adults %>%
  mutate(Species_Name = case_when(
    Species_Name == "Holhymenia_tibialis" ~ "Holhymenia_tibialis / clavigera",
    TRUE ~ Species_Name
  )) %>% 
  filter(Antena_Adulto != 0, 
         !Species_Name %in% c("Jadera_haematoloma", "Harmostes_serratus",
                              "Oedancala_sp.", "Oncopeltus_fasciatus",
                              "Largus_sp.", "Dysdercus_mimus",
                              "Dysdercus_suturellus", "Halyomorpha_halys"),
         Species_Name %in% dados.ninfa$Species_Name)

# Now we prune the tree to the data we have on Coreidae adults genera

t_adults_gen = keep.tip(t1, adults_data2$Species_Name)

# adults matrix

adults.gen.sim = setNames(as.factor(adults_data2$Soma), adults_data2$Species_Name)
adults.gen.sim

# Mk model
#fit.ER.equal.ag = fitMk(ladderize(t_adults_gen), adults.gen.sim, model = "ER",
#           pi = 'equal')


#fit.ER.fitz.ag = fitMk(ladderize(t_adults_gen), adults.gen.sim , model = "ER",
#                 pi = 'fitzjohn')

fit.ER.equal.ag
fit.ER.fitz.ag

plot(fit.ER.equal.ag)

save(fit.ER.equal.ag, file = './evo.models/er-equal-adults-genus.RData')
save(fit.ER.fitz.ag, file = './evo.models/er-fitz-adults-genus.RData')

load( file = './evo.models/er-equal-adults-genus.RData')
load( file = './evo.models/er-fitz-adults-genus.RData')

#fit.ARD.equal.ag = fitMk(ladderize(t_adults_gen), adults.gen.sim, model = "ARD",
#                        pi = 'equal')

#fit.ARD.fitz.ag = fitMk(ladderize(t_adults_gen), adults.gen.sim, model = "ARD",
#                       pi = 'fitzjohn')

fit.ARD.equal.ag
fit.ARD.fitz.ag

save(fit.ARD.equal.ag, file = './evo.models/ARD-equal-adults-genus.RData')
save(fit.ARD.fitz.ag, file = './evo.models/ARD-fitz-adults-genus.RData')

load( file = './evo.models/ARD-equal-adults-genus.RData')
load( file = './evo.models/ARD-fitz-adults-genus.RData')

# Hidden rates model, HRM

#fitHRM.er.equal.ag = fitHRM(ladderize(t_adults_gen), adults.gen.sim, model = "ER",
#                           pi = 'equal')

#fitHRM.er.fitz.ag = fitHRM(ladderize(t_adults_gen), adults.gen.sim, model = "ER",
#                          pi = 'fitzjohn')

save(fitHRM.er.equal.ag, file = './evo.models/hrm-er-equal-adults-genus.RData')
save(fitHRM.er.fitz.ag, file = './evo.models/hrm-er-fitz-adults-genus.RData')

load(file = './evo.models/hrm-er-equal-adults-genus.RData')
load(file = './evo.models/hrm-er-fitz-adults-genus.RData')

plot(fitHRM.er.equal.ag)
fitHRM.er.equal.ag

plot(fitHRM.er.fitz.ag)
fitHRM.er.fitz.ag

#fitHRM.ard.equal.ag = fitHRM(ladderize(t_adults_gen), adults.gen.sim, model = "ARD",
#                            pi = 'equal')


#fitHRM.ard.fitz.ag = fitHRM(ladderize(t_adults_gen), adults.gen.sim, model = "ARD",
#                           pi = 'fitzjohn')

save(fitHRM.ard.equal.ag, file = './evo.models/hrm-equal-adults-genus.RData')
save(fitHRM.ard.fitz.ag, file = './evo.models/hrm-fitz-adults-genus.RData')


load(file = './evo.models/hrm-equal-adults-genus.RData')
load(file = './evo.models/hrm-fitz-adults-genus.RData')


# Comparing models
fit_models_ag = anova(fit.ER.equal.ag, fit.ARD.equal.ag,
                      fitHRM.er.equal.ag, fitHRM.ard.equal.ag,
                      fit.ER.fitz.ag, fit.ARD.fitz.ag,
                      fitHRM.er.fitz.ag, fitHRM.ard.fitz.ag)

# It was chosen the ER fitz model once it ad the lowest AIC value abd
# degree of freedom of 1

fit.ER.fitz.ag
plot(fit.ER.fitz.ag)


#### PLOTTING NYMPHS AND ADULTS GENERA TREE TOGETHER ####

#Nymph 
anc.nymph = ancr(fitHRM.er.fitz.nymph, hide.hidden = F)

#Adults
anc.ag = ancr(fit.ER.fitz.ag)

tree_plot = fit.ER.fitz.ag$tree

anc.nymph$tree = tree_plot
anc.ag$tree = tree_plot

states = colnames(anc.nymph$ace)
cols = setNames(
  c("#0E3768","blue", "#94129F", "#ea22fb",  
    "#DE627A", "pink", "#ce4326", "#ff8b5f"),
  states
)

states2 = colnames(anc.ag$ace)
cols2 = setNames(
  c("#0E3768", "#94129F",
    "#DE627A",  "#ce4326"),
  states2
)
cols2

tip_state = as.character(nymphs.gen.sim[anc.nymph$tree$tip.label])
states2

tip.pies = sapply(
  states2,
  function(z) as.numeric(tip_state == z)
)

tip.pies = as.matrix(tip.pies)
rownames(tip.pies) = anc.nymph$tree$tip.label
colnames(tip.pies) = states2

getwd()
png(file = './figures/ancestral-facing.png', height = 400, width = 350,
    units = "mm", res = 300)
pdf(file = './figures/antennnomere_facing.pdf',
    h = 12, w = 12)

layout(matrix(1:3, nrow = 1), widths = c(0.44, 0.12, 0.44))

plot(
  anc.nymph,
  args.plotTree = list(
    ftype = "off",
    lwd = 1.5,
    ylim = c(0, Ntip(tree_plot) + 1),
    mar = c(1, 0.5, 3, 0.5)
  ),
  args.nodelabels = list(
    piecol = cols,
    cex = 0.3
  )
)
par(fg = 'transparent')
tiplabels(
  pie = tip.pies,
  piecol = cols2,
  cex = 0.45,
  border = "transparent"
)
par(fg = 'black')
mtext("nymph", side = 3, line = 1, adj = 0)

## Middle: taxon labels, plotted once
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

## Tips are nodes 1:Ntip(tree); arrange bottom to top
tip_order <- order(pp$yy[seq_len(Ntip(tree_plot))])

tip_labels <- tree_plot$tip.label[tip_order]
tip_y      <- pp$yy[tip_order]
y_lim      <- pp$y.lim

## Middle panel: use the plotted y coordinates, NOT 1:Ntip()
plot.new()
plot.window(xlim = c(-1, 1), ylim = y_lim)

text(
  x = 0,
  y = tip_y,
  labels = gsub("_", " ", tip_labels),
  cex = 1,
  font = 3
)


plot(
  anc.ag,
  args.plotTree = list(
    ftype = "off",
    direction = "leftwards",
    lwd = 1.5,
    ylim = c(0, Ntip(tree_plot) + 1),
    mar = c(1, 0.5, 3, 0.5)
  ),
  args.nodelabels = list(
    piecol = cols2,
    cex = 0.3
  ),
  args.tiplabels = list(
    cex = 0.45
  )
)
mtext("adults", side = 3, line = 1, adj = 1)
dev.off()

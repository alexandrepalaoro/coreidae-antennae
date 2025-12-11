library(ape)
library(phytools)
library(geiger)
library(diversitree)

#### LOADING PHYLOGENY AND DATA FILE ####

tt = read.tree(file = './phylo/50p.treefile')

tips = read.csv(file = './data/23-10-25_antena_coreidae_summarised_2.csv', h =T, sep = ';')

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


#### DATING THE PHYLOGENY #####

# There are several methods for dating the tree in the chronos function 
# We intended to perform analyses using at least three types to ensure results 
# were robust regardless of branch lengths. However, we realized that changes 
# in antennal state are relativelly rare in the phylogeny. So, we left this
# here to show that we intended to do it, but the patterns found were robust
# to not need such a sensitivity analyses.

## This is the data sheet containing the age of the nodes according to the fossil
## records.
calib.t2 = read.csv(file = './data/29-04-25_calib_points_PBD.csv', h = T, sep = ';')

# Now, the chronos analyses

### NOTE: I COMMENTED THE FUNCTIONS BECAUSE I SAVED THE RESULT.
### THUS, YOU CAN LOAD THE FILES INSTEAD OF RUNNING THE ANALYSES YOURSELF.
### IF YOU WANT TO RUN THE FUNCTION, JUST UNCOMMENT THE FUNCTIONS.

#corr.t2 = chronos(t1, lambda = 1, model = 'correlated', 
#                 calibration = calib.t2)

#save(corr.t2, file = './corr-t2-antena.RData')
load(file = "./corr-t2-antena.RData")

plot(corr.t2)

#clock.t2 = chronos(t1, lambda = 0.565, model = 'clock', 
#                   calibration = calib.t2)
#save(clock.t2, file = './clock-t2-antena.RData')
#load(file = "./clock-t2-antena.RData")

#plot(clock.t2)


#relax.t2 = chronos(t1, lambda = 1, model = 'relaxed', 
#                  calibration = calib.t2)
#save(relax.t2, file = './relax-t2-antena.RData')
load(file = "./relax-t2-antena.RData")

plot(relax.t2)

#### REMOVING ZEROES FROM DATASET AND DROPPING UNNECESSARY TIPS ####

tips$Antena

#antena = tips


ant = read.csv(file = './data/23-10-25_Coreidae_antena.csv', h =T, sep = ';')
rownames(ant) = ant$Species_Name

#antena = tips[tips$Antena!="No_data",]
#antena$Antena

# Since our goal is analyzing trait evolution in Coreidae,
# we removed the non-Coreidae species in the tree.

antena = ant[ant$Species_Name != "Jadera_haematoloma",]
antena = ant[ant$Species_Name != "Harmostes_serratus",]
antena = ant[ant$Species_Name != "Oedancala_sp.",]
antena = ant[ant$Species_Name != "Oncopeltus_fasciatus",]
antena = ant[ant$Species_Name != "Largus_sp.",]
antena = ant[ant$Species_Name != "Dysdercus_mimus",]
antena = ant[ant$Species_Name != "Dysdercus_suturellus",]
antena = ant[ant$Species_Name != "Halyomorpha_halys",]

# Now we prune the tree to the data we have on Coreidae.

#clock.t2.ant = keep.tip(clock.t2, tip = antena$Species_Name)

corr.t2.ant = keep.tip(corr.t2, tip = antena$Species_Name)

relax.t2.ant = keep.tip(relax.t2, tip = antena$Species_Name)

#name.check(clock.t2.ant,antena)
name.check(corr.t2.ant,antena)
name.check(relax.t2.ant,antena)


#### THE MK MODEL: ASSESSING CHANGES IN ANTENNAL TRAIT STATE ####

# Before running the analyses, we need to create a matrix so phytools
# can read our data.


ant.sim = setNames(as.factor(antena$Antena), rownames(antena))

ant.sim

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

#fit.ER.equal = fitMk(ladderize(corr.t2.ant), ant.sim , model = "ER",
#                     pi = 'equal')


#fit.ER.fitz = fitMk(ladderize(corr.t2.ant), ant.sim , model = "ER",
#                    pi = 'fitzjohn')

#plot(fit.ER.fitz)

#save(fit.ER.equal, file = './evo.models/er-equal-antena.RData')
#save(fit.ER.fitz, file = './evo.models/er-fitz-antena.RData')

load( file = './evo.models/er-equal-antena.RData')
load( file = './evo.models/er-fitz-antena.RData')

plot(fit.ER.equal)
fit.ER.equal

plot(fit.ER.fitz)
fit.ER.fitz

# All rates differ model #

#fit.ARD.equal = fitMk(ladderize(corr.t2.ant), ant.sim, model = "ARD",
#                      pi = 'equal')

#fit.ARD.fitz = fitMk(ladderize(corr.t2.ant), ant.sim, model = "ARD",
#                     pi = 'fitzjohn')

#save(fit.ARD.equal, file = './evo.models/ard-equal-antena.RData')
#save(fit.ARD.fitz, file = './evo.models/ard-fitz-antena.RData')

load(file = './evo.models/ard-equal-antena.RData')
load(file = './evo.models/ard-fitz-antena.RData')

plot(fit.ARD.equal)
fit.ARD.equal

plot(fit.ARD.fitz)
fit.ARD.fitz

# We are also using the HRM models for thoroughness. Maybe there is an
# underlying mechanism for changing antennal state we are not capturing. 

#fitHRM.er.equal = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ER",
#                         pi = 'equal')

#fitHRM.er.fitz = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ER",
#                        pi = 'fitzjohn')

#save(fitHRM.er.equal, file = './evo.models/hrm-er-equal-antena.RData')
#save(fitHRM.er.fitz, file = './evo.models/hrm-er-fitz-antena.RData')

load(file = './evo.models/hrm-er-equal-antena.RData')
load(file = './evo.models/hrm-er-fitz-antena.RData')

plot(fitHRM.er.equal)
fitHRM.er.equal

plot(fitHRM.er.fitz)
fitHRM.er.fitz

# All rates differ HRM #

#fitHRM.equal = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ARD",
#                      pi = 'equal')


#fitHRM.fitz = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ARD",
#                     pi = 'fitzjohn')

#save(fitHRM.equal, file = './evo.models/hrm-equal-antena.RData')
#save(fitHRM.fitz, file = './evo.models/hrm-fitz-antena.RData')


load(file = './evo.models/hrm-equal-antena.RData')
load(file = './evo.models/hrm-fitz-antena.RData')


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
  
anc.ARD = ancr(fit.ARD.equal)

fit.ARD.equal

cols = setNames(c("#DE627A","#0E3768"),
               levels(ant.sim))


pdf(file = './figures/ancestral-tree-fan.pdf', h = 24, w = 20)
plot(anc.ARD,  
     args.plotTree = list(lwd = 2, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          type = 'fan',
                          cex = 0.5, offset = 3
                          ),
     args.nodelabels = list(cex = 0.4, piecol = cols),
     args.tiplabels = list(cex = 0.3)
     
)
dev.off()


pdf(file = './figures/ancestral-tree-phylogram.pdf', h = 26, w = 16)
plot(anc.ARD,  
     args.plotTree = list(lwd = 2, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          type = 'phylogram',
                          cex = 0.5, offset = 0.6),
     args.nodelabels = list(cex = 0.25, piecol = cols),
     args.tiplabels = list(cex = 0.2)
     
)
dev.off()



ard_sim = make.simmap( tree = ladderize(corr.t2.ant), x= ant.sim , model = "ARD", pi = 'equal',
                     nsim = 1000)

sum_sim = describe.simmap(ard_sim)
sum_sim

plot(density(ard_sim))


#### PLOT FULL TREE WITH MISSING DATA ####

# This is used to plot the entire tree from Forthman et al. (2024)
# and showing which species we have data on (and in this case, 
# what is the antennal steate), and the species we do not have data
# on. 

full = read.csv(file = './data/29-04-25_Coreidae_antena_uncertain_2.csv', h =T, sep = ';')
rownames(full) = full$Species_Name

corr.t2.full = keep.tip(corr.t2, tip = full$Species_Name)

name.check(corr.t2.full,full)

antena.full = data.frame(Antena = as.factor(full$Antenna))
rownames(antena.full) = rownames(full)

name.check(corr.t2.full,antena.full)


colors<-list(
  setNames(c("#DE627A","white", "#0E3768"),c("Expansion","No_data","Straight")))

pdf(file = './figures/full-tree.pdf', h = 10, w = 18)
t = plotFanTree.wTraits(corr.t2.full, antena.full,
                        colors = colors,
                        spacer = 0,
                        fsize = 0.7, part = 0.5)
legend(x=144.4429,y=453.3577,
       names(colors[[1]]),border = 'black',
       fill=colors[[1]],
       title="Antennal morphology",bty="n",xjust=0.5,yjust=0.5)
dev.off()  



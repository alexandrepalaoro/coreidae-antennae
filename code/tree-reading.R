library(ape)
library(phytools)
library(geiger)
library(diversitree)

tt = read.tree(file = './phylo/50p.treefile')

#x = tt$tip.label

#write.csv(tt$tip.label, file = './data/tips2.csv' )

tips = read.csv(file = './data/antena_coreidae_3.csv', h =T, sep = ';')

tips = read.csv(file = './data/antena_coreidae_summarised_2.csv', h =T, sep = ';')
tips = read.csv(file = './data/antena_coreidae_summarised.csv', h =T, sep = ';')

rownames(tips) = tips$Species_Name

tt$tip.label = tips$Species_Name

plot(tt)

name.check(tt, tips)


#tips.ninfa = read.csv(file = './data/antena_ninfa_coreidae.csv', h = T)
#rownames(tips.ninfa) = tips.ninfa$Species_name

#name.check(tt,tips.ninfa)

#t1 = root(tt, outgroup =  "Jadera_haematoloma", resolve.root = F)

t1 = root(tt, node = 308, resolve.root = T)

plot.phylo(ladderize(t1))
plot.phylo(t1)

#t2 = drop.tip(t1, "Jadera_haematoloma")

#plot(t2)


#### DATING THE PHYLOGENY

#lad.t1 = ladderize(t1)

#plot(lad.t1)
#nodelabels()

pdf(file = './tree-for-calib.pdf', h = 40, w = 35)
plot.phylo(t1)
nodelabels()
dev.off()

calib.t1 = read.csv(file = './data/calib_points_PBD.csv', h = T, sep = ';')

#teste = chronos(t1, lambda = 1, model = 'relaxed', 
#                calibration = calib.t1)
#teste = chronos(t1, model = 'clock')

#plot(teste)


calib.t2 = read.csv(file = './data/calibration_points_non-ladder.csv', h = T, sep = ';')

corr.t2 = chronos(t1, lambda = 1, model = 'correlated', 
                calibration = calib.t1)
#save(corr.t2, file = './corr-t2-antena.RData')
load(file = "./corr-t2-antena.RData")

plot(corr.t2)


clock.t2 = chronos(t1, lambda = 1, model = 'clock', 
                   calibration = calib.t1)
#save(clock.t2, file = './clock-t2-antena.RData')
load(file = "./clock-t2-antena.RData")

plot(clock.t2)


relax.t2 = chronos(t1, lambda = 1, model = 'relaxed', 
                   calibration = calib.t1)
#save(relax.t2, file = './relax-t2-antena.RData')
load(file = "./relax-t2-antena.RData")

plot(relax.t2)

#### ADULTS ####
#### REMOVING ZEROES FROM DATASET AND DROPPING UNNECESSARY TIPS ####

tips$Antena

#antena = tips

antena = tips[tips$Antena!="0",]
antena$Antena
antena = antena[antena$Species_Name != "Jadera_haematoloma",]
antena = antena[antena$Species_Name != "Harmostes_serratus",]

antena = antena[antena$Species_Name != "Oedancala_sp.",]
antena = antena[antena$Species_Name != "Oncopeltus_fasciatus",]
antena = antena[antena$Species_Name != "Largus_sp.",]
antena = antena[antena$Species_Name != "Dysdercus_mimus",]
antena = antena[antena$Species_Name != "Dysdercus_suturellus",]
antena = antena[antena$Species_Name != "Halyomorpha_halys",]



clock.t2.ant = keep.tip(clock.t2, tip = antena$Species_Name)

corr.t2.ant = keep.tip(corr.t2, tip = antena$Species_Name)

relax.t2.ant = keep.tip(relax.t2, tip = antena$Species_Name)

name.check(clock.t2.ant,antena)
name.check(corr.t2.ant,antena)
name.check(relax.t2.ant,antena)


plot(corr.t2.ant)

#### CHANGES IN STATES ####

ant.sim = setNames(as.factor(antena$Antena), rownames(antena))

#trees = make.simmap(ladderize(clock.t2.ant),ant.sim, model = "ER",nsim = 100)
#trees1 = make.simmap(ladderize(clock.t2.ant), ant.sim, model = "ARD", nsim = 1000)

#obj2 = summary(trees1)

fit.ER.null = fitMk(ladderize(corr.t2.ant), ant.sim , model = "ER", 
                    pi = c(0.5,0.5,0,0))
plot(fit.ER.null)

fit.ER.equal = fitMk(ladderize(corr.t2.ant), ant.sim , model = "ER",
               pi = 'equal')
plot(fit.ER.equal)

fit.ER.fitz = fitMk(ladderize(corr.t2.ant), ant.sim , model = "ER",
                     pi = 'fitzjohn')
plot(fit.ER.fitz)

save(fit.ER.null, file = './evo.models/er-null-antena.RData')
save(fit.ER.equal, file = './evo.models/er-equal-antena.RData')
save(fit.ER.fitz, file = './evo.models/er-fitz-antena.RData')

load( file = './evo.models/er-null-antena.RData')
load(file = './evo.models/er-equal-antena.RData')
load( file = './evo.models/er-fitz-antena.RData')



fit.ARD.null = fitMk(ladderize(corr.t2.ant), ant.sim, model = "ARD",
                     pi = c(0.5,0.5,0,0))
plot(fit.ARD.null)

fit.ARD.equal = fitMk(ladderize(corr.t2.ant), ant.sim, model = "ARD",
                     pi = 'equal')
plot(fit.ARD.equal)

fit.ARD.fitz = fitMk(ladderize(corr.t2.ant), ant.sim, model = "ARD",
                      pi = 'fitzjohn')
plot(fit.ARD.fitz)

save(fit.ARD.null, file = './evo.models/ard-null-antena.RData')
save(fit.ARD.equal, file = './evo.models/ard-equal-antena.RData')
save(fit.ARD.fitz, file = './evo.models/ard-fitz-antena.RData')


load(file = './evo.models/ard-null-antena.RData')
load(file = './evo.models/ard-equal-antena.RData')
load(file = './evo.models/ard-fitz-antena.RData')

fitHRM.er.null = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ER",
                        pi = c(0.5,0,0.5,0,0,0,0,0))
plot(fitHRM.er.null)


fitHRM.er.equal = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ER",
                        pi = 'equal')
plot(fitHRM.er.equal)


fitHRM.er.fitz = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ER",
                         pi = 'fitzjohn')
plot(fitHRM.er.fitz)

save(fitHRM.er.null, file = './evo.models/hrm-er-null-antena.RData')
save(fitHRM.er.equal, file = './evo.models/hrm-er-equal-antena.RData')
save(fitHRM.er.fitz, file = './evo.models/hrm-er-fitz-antena.RData')

load(file = './evo.models/hrm-er-null-antena.RData')
load(file = './evo.models/hrm-er-equal-antena.RData')
load(file = './evo.models/hrm-er-fitz-antena.RData')


fitHRM.null = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ARD",
                     pi = c(0.5,0,0.5,0,0,0,0,0))
plot(fitHRM.null)


fitHRM.equal = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ARD",
                     pi = 'equal')
plot(fitHRM.equal)


fitHRM.fitz = fitHRM(ladderize(corr.t2.ant), ant.sim , model = "ARD",
                      pi = 'fitzjohn')
plot(fitHRM.fitz)

save(fitHRM.null, file = './evo.models/hrm-null-antena.RData')
save(fitHRM.equal, file = './evo.models/hrm-equal-antena.RData')
save(fitHRM.fitz, file = './evo.models/hrm-fitz-antena.RData')

load(file = './evo.models/hrm-null-antena.RData')
load(file = './evo.models/hrm-equal-antena.RData')
load(file = './evo.models/hrm-fitz-antena.RData')


AIC(fit.ER.null, fit.ARD.null,
    fitHRM.er.null, fitHRM.null,
    fit.ER.equal, fit.ARD.equal,
    fitHRM.er.equal, fitHRM.equal,
    fit.ER.fitz, fit.ARD.fitz,
    fitHRM.er.fitz, fitHRM.fitz
    )

AIC(fit.ER.equal, fit.ARD.null,
    fitHRM.er.null, fitHRM.null)

AIC(fit.ARD.corr.null) - AIC(fitHRM.er.null)

plot(fit.ER.fitz)
print(fit.ER.fitz)
anc.ARD = ancr(fit.ER.fitz)

cols<-setNames(c("gold1","olivedrab3","darkgreen"),
               levels(ant.sim))

pdf(file = './ancestral-tree-fan.pdf', h = 28, w = 20)
plot(anc.ARD,  
     args.plotTree = list(lwd = 2, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          type = 'fan',
                          cex = 0.5, offset = 3),
     args.nodelabels = list(cex = 0.4),
     args.tiplabels = list(cex = 0.3)
     
     )
dev.off()


pdf(file = './ancestral-tree-phylogram.pdf', h = 24, w = 16)
plot(anc.ARD,  
     args.plotTree = list(lwd = 2, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          type = 'phylogram',
                          cex = 0.5, offset = 0.6),
     args.nodelabels = list(cex = 0.4),
     args.tiplabels = list(cex = 0.2)
     
)
dev.off()



anc.ARD.equal = ancr(fit.ARD.equal)

pdf(file = './ancestral-tree-equal.pdf', h = 28, w = 20)
plot(anc.ARD.equal, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()



plot(fit.ARD.corr.null)
anc.ard.corr = ancr(fit.ARD.corr.null)

pdf(file = './ancestral-tree-corr.pdf', h = 28, w = 20)
plot(anc.ard.corr, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()



#### NYMPHS ####
#### REMOVING ZEROES FROM DATASET AND DROPPING UNNECESSARY TIPS ####

nymph = read.csv(file = './data/antena_ninfa.csv', h = T, sep = ';')

rownames(nymph) = nymph$Species_Name

nymph$Antena

ant.nymph = nymph[nymph$Antena!="0",]
ant.nymph$Antena
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Jadera_haematoloma",]
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Harmostes_serratus",]

ant.nymph = ant.nymph[ant.nymph$Species_Name != "Oedancala_sp.",]
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Oncopeltus_fasciatus",]
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Largus_sp.",]
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Dysdercus_mimus",]
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Dysdercus_suturellus",]
ant.nymph = ant.nymph[ant.nymph$Species_Name != "Halyomorpha_halys",]

clock.t2.nymph = keep.tip(clock.t2, tip = ant.nymph$Species_Name)

corr.t2.nymph = keep.tip(corr.t2, tip = ant.nymph$Species_Name)

relax.t2.nymph = keep.tip(relax.t2, tip = ant.nymph$Species_Name)

name.check(clock.t2.nymph,ant.nymph)
name.check(corr.t2.nymph,ant.nymph)
name.check(relax.t2.nymph,ant.nymph)


#### MODELS OF EVOLUTION ####.

nymph.sim = setNames(as.factor(ant.nymph$Antena), rownames(ant.nymph))


nymph.ER.null = fitMk(ladderize(corr.t2.nymph), nymph.sim , model = "ER", 
                    pi = c(0.5,0.5,0,0))
plot(nymph.ER.null)

nymph.ER.equal = fitMk(ladderize(corr.t2.nymph), nymph.sim , model = "ER",
                     pi = 'equal')
plot(nymph.ER.equal)

nymph.ER.fitz = fitMk(ladderize(corr.t2.nymph), nymph.sim , model = "ER",
                    pi = 'fitzjohn')
plot(nymph.ER.fitz)

save(nymph.ER.null, file = './evo.models/er-null-nymph.RData')
save(nymph.ER.equal, file = './evo.models/er-equal-nymph.RData')
save(nymph.ER.fitz, file = './evo.models/er-fitz-nymp.RData')

load( file = './evo.models/er-null-nymph.RData')
load( file = './evo.models/er-equal-nymph.RData')
load( file = './evo.models/er-fitz-nymp.RData')


nymph.ARD.null = fitMk(ladderize(corr.t2.nymph), nymph.sim, model = "ARD",
                     pi = c(0.5,0.5,0,0))
plot(nymph.ARD.null)

nymph.ARD.equal = fitMk(ladderize(corr.t2.nymph), nymph.sim, model = "ARD",
                      pi = 'equal')
plot(nymph.ARD.equal)

nymph.ARD.fitz = fitMk(ladderize(corr.t2.nymph), nymph.sim, model = "ARD",
                     pi = 'fitzjohn')
plot(nymph.ARD.fitz)

save(nymph.ARD.null, file = './evo.models/ard-null-nymph.RData')
save(nymph.ARD.equal, file = './evo.models/ard-equal-nymph.RData')
save(nymph.ARD.fitz, file = './evo.models/ard-fitz-nymph.RData')


load(file = './evo.models/ard-null-nymph.RData')
load(file = './evo.models/ard-equal-nymph.RData')
load(file = './evo.models/ard-fitz-nymph.RData')

nymph.fitHRM.er.null = fitHRM(ladderize(corr.t2.nymph), nymph.sim , model = "ER",
                        pi = c(0.5,0,0.5,0,0,0,0,0))
plot(nymph.fitHRM.er.null)


nymph.fitHRM.er.equal = fitHRM(ladderize(corr.t2.nymph), nymph.sim , model = "ER",
                         pi = 'equal')
plot(nymph.fitHRM.er.equal)


nymph.fitHRM.er.fitz = fitHRM(ladderize(corr.t2.nymph), nymph.sim , model = "ER",
                        pi = 'fitzjohn')
plot(nymph.fitHRM.er.fitz)

save(nymph.fitHRM.er.null, file = './evo.models/hrm-er-null-nymph.RData')
save(nymph.fitHRM.er.equal, file = './evo.models/hrm-er-equal-nymph.RData')
save(nymph.fitHRM.er.fitz, file = './evo.models/hrm-er-fitz-nymph.RData')

load(file = './evo.models/hrm-er-null-nymph.RData')
load(file = './evo.models/hrm-er-equal-nymph.RData')
load(file = './evo.models/hrm-er-fitz-nymph.RData')


nymph.fitHRM.null = fitHRM(ladderize(corr.t2.nymph), nymph.sim , model = "ARD",
                     pi = c(0.5,0,0.5,0,0,0,0,0))
plot(nymph.fitHRM.null)


nymph.fitHRM.equal = fitHRM(ladderize(corr.t2.nymph), nymph.sim , model = "ARD",
                      pi = 'equal')
plot(nymph.fitHRM.equal)


nymph.fitHRM.fitz = fitHRM(ladderize(corr.t2.nymph), nymph.sim , model = "ARD",
                     pi = 'fitzjohn')
plot(nymph.fitHRM.fitz)

save(nymph.fitHRM.null, file = './evo.models/hrm-null-nymph.RData')
save(nymph.fitHRM.equal, file = './evo.models/hrm-equal-nymph.RData')
save(nymph.fitHRM.fitz, file = './evo.models/hrm-fitz-nymph.RData')

load(file = './evo.models/hrm-null-nymph.RData')
load(file = './evo.models/hrm-equal-nymph.RData')
load(file = './evo.models/hrm-fitz-nymph.RData')


AIC(nymph.ER.null, nymph.ARD.null,
    nymph.fitHRM.er.null, nymph.fitHRM.null,
    nymph.ER.equal, nymph.ARD.equal,
    nymph.fitHRM.er.equal, nymph.fitHRM.equal,
    nymph.ER.fitz, nymph.ARD.fitz,
    nymph.fitHRM.er.fitz, nymph.fitHRM.fitz
)


anc.nymph.ER = ancr(nymph.fitHRM.er.fitz)

plot(nymph.fitHRM.er.fitz)


pdf(file = './ancestral-nymph-tree.pdf', h = 28, w = 20)
plot(anc.nymph.ER, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()



### ADULT GENERA ONLY ###

adult_genera = read.csv(file = './data/Coreidae_antena - Adultos_genero.csv', h = T, sep = ';')

adult_anl = adult_genera[adult_genera$Antena!='0',] 

rownames(adult_anl) = adult_anl$Species_Name

genera_adult_clock = keep.tip(clock.t2, tip = adult_anl$Species_Name)
name.check(genera_adult_clock,adult_anl)

ant.adult.sim = setNames(as.factor(adult_anl$Antena), rownames(adult_anl))

fit.ER.null = fitMk(ladderize(genera_adult_clock), ant.adult.sim , model = "ER", 
                    pi = c(0.5,0.5,0,0))
plot(fit.ER.null)

fit.ER.equal = fitMk(ladderize(genera_adult_clock), ant.adult.sim , model = "ER",
                     pi = 'equal')
plot(fit.ER.equal)

fit.ER.fitz = fitMk(ladderize(genera_adult_clock), ant.adult.sim , model = "ER",
                    pi = 'fitzjohn')
plot(fit.ER.fitz)

save(fit.ER.null, file = './evo.models/adult-er-null-antena.RData')
save(fit.ER.equal, file = './evo.models/adult-er-equal-antena.RData')
save(fit.ER.fitz, file = './evo.models/adult-er-fitz-antena.RData')

load( file = './evo.models/adult-er-null-antena.RData')
load(file = './evo.models/adult-er-equal-antena.RData')
load( file = './evo.models/adult-er-fitz-antena.RData')



fit.ARD.null = fitMk(ladderize(genera_adult_clock), ant.adult.sim, model = "ARD",
                     pi = c(0.5,0.5,0,0))
plot(fit.ARD.null)

fit.ARD.equal = fitMk(ladderize(genera_adult_clock), ant.adult.sim, model = "ARD",
                      pi = 'equal')
plot(fit.ARD.equal)

fit.ARD.fitz = fitMk(ladderize(genera_adult_clock), ant.adult.sim, model = "ARD",
                     pi = 'fitzjohn')
plot(fit.ARD.fitz)

save(fit.ARD.null, file = './evo.models/adult-ard-null-antena.RData')
save(fit.ARD.equal, file = './evo.models/adult-ard-equal-antena.RData')
save(fit.ARD.fitz, file = './evo.models/adult-ard-fitz-antena.RData')


load(file = './evo.models/adult-ard-null-antena.RData')
load(file = './evo.models/adult-ard-equal-antena.RData')
load(file = './evo.models/adult-ard-fitz-antena.RData')

fitHRM.er.null = fitHRM(ladderize(genera_adult_clock), ant.adult.sim , model = "ER",
                        pi = c(0.5,0,0.5,0,0,0,0,0))
plot(fitHRM.er.null)


fitHRM.er.equal = fitHRM(ladderize(genera_adult_clock), ant.adult.sim , model = "ER",
                         pi = 'equal')
plot(fitHRM.er.equal)


fitHRM.er.fitz = fitHRM(ladderize(genera_adult_clock), ant.adult.sim , model = "ER",
                        pi = 'fitzjohn')
plot(fitHRM.er.fitz)

save(fitHRM.er.null, file = './evo.models/adult-hrm-er-null-antena.RData')
save(fitHRM.er.equal, file = './evo.models/adult-hrm-er-equal-antena.RData')
save(fitHRM.er.fitz, file = './evo.models/adult-hrm-er-fitz-antena.RData')

load(file = './evo.models/adult-hrm-er-null-antena.RData')
load(file = './evo.models/adult-hrm-er-equal-antena.RData')
load(file = './evo.models/adult-hrm-er-fitz-antena.RData')


fitHRM.null = fitHRM(ladderize(genera_adult_clock), ant.adult.sim , model = "ARD",
                     pi = c(0.5,0,0.5,0,0,0,0,0))
plot(fitHRM.null)


fitHRM.equal = fitHRM(ladderize(genera_adult_clock), ant.adult.sim , model = "ARD",
                      pi = 'equal')
plot(fitHRM.equal)


fitHRM.fitz = fitHRM(ladderize(genera_adult_clock), ant.adult.sim , model = "ARD",
                     pi = 'fitzjohn')
plot(fitHRM.fitz)

save(fitHRM.null, file = './evo.models/adult-hrm-null-antena.RData')
save(fitHRM.equal, file = './evo.models/adult-hrm-equal-antena.RData')
save(fitHRM.fitz, file = './evo.models/adult-hrm-fitz-antena.RData')

load(file = './evo.models/adult-hrm-null-antena.RData')
load(file = './evo.models/adult-hrm-equal-antena.RData')
load(file = './evo.models/adult-hrm-fitz-antena.RData')


AIC(fit.ER.null, fit.ARD.null,
    fitHRM.er.null, fitHRM.null,
    fit.ER.equal, fit.ARD.equal,
    fitHRM.er.equal, fitHRM.equal,
    fit.ER.fitz, fit.ARD.fitz,
    fitHRM.er.fitz, fitHRM.fitz
)

AIC(fit.ER.equal, fit.ARD.null,
    fitHRM.er.null, fitHRM.null)

AIC(fit.ARD.corr.null) - AIC(fitHRM.er.null)

plot(fit.ARD.equal)
print(fit.ARD.equal)
anc.ARD = ancr(fit.ARD.equal)

pdf(file = './ancestral-tree-adult-genera.pdf', h = 28, w = 20)
plot(anc.ARD, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()

anc.ARD.null = ancr(fit.ARD.null)

pdf(file = './ancestral-tree-null.pdf', h = 28, w = 20)
plot(anc.ARD.null, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()



plot(fit.ARD.corr.null)
anc.ard.corr = ancr(fit.ARD.corr.null)

pdf(file = './ancestral-tree-corr.pdf', h = 28, w = 20)
plot(anc.ard.corr, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()


#### NYMPH GENERA ONLY ####

nymph_genera = read.csv(file = './data/Coreidae_antena - Ninfa_genero.csv', h = T,
                        sep = ';')

nymph_anl = nymph_genera[nymph_genera$Antena!='0',] 

rownames(nymph_anl) = nymph_anl$Species_Name

genera_nymph_clock = keep.tip(clock.t2, tip = nymph_anl$Species_Name)

name.check(genera_nymph_clock,nymph_anl)

ant.nymph.sim = setNames(as.factor(nymph_anl$Antena), rownames(nymph_anl))

fit.ER.null = fitMk(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ER", 
                    pi = c(0.5,0.5,0,0))
plot(fit.ER.null)

fit.ER.equal = fitMk(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ER",
                     pi = 'equal')
plot(fit.ER.equal)

fit.ER.fitz = fitMk(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ER",
                    pi = 'fitzjohn')
plot(fit.ER.fitz)

save(fit.ER.null, file = './evo.models/nymph-er-null-antena.RData')
save(fit.ER.equal, file = './evo.models/nymph-er-equal-antena.RData')
save(fit.ER.fitz, file = './evo.models/nymph-er-fitz-antena.RData')

load( file = './evo.models/nymph-er-null-antena.RData')
load(file = './evo.models/nymph-er-equal-antena.RData')
load( file = './evo.models/nymph-er-fitz-antena.RData')



fit.ARD.null = fitMk(ladderize(genera_nymph_clock), ant.nymph.sim, model = "ARD",
                     pi = c(0.5,0.5,0,0))
plot(fit.ARD.null)

fit.ARD.equal = fitMk(ladderize(genera_nymph_clock), ant.nymph.sim, model = "ARD",
                      pi = 'equal')
plot(fit.ARD.equal)

fit.ARD.fitz = fitMk(ladderize(genera_nymph_clock), ant.nymph.sim, model = "ARD",
                     pi = 'fitzjohn')
plot(fit.ARD.fitz)

save(fit.ARD.null, file = './evo.models/nymph-ard-null-antena.RData')
save(fit.ARD.equal, file = './evo.models/nymph-ard-equal-antena.RData')
save(fit.ARD.fitz, file = './evo.models/nymph-ard-fitz-antena.RData')


load(file = './evo.models/nymph-ard-null-antena.RData')
load(file = './evo.models/nymph-ard-equal-antena.RData')
load(file = './evo.models/nymph-ard-fitz-antena.RData')

fitHRM.er.null = fitHRM(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ER",
                        pi = c(0.5,0,0.5,0,0,0,0,0))
plot(fitHRM.er.null)


fitHRM.er.equal = fitHRM(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ER",
                         pi = 'equal')
plot(fitHRM.er.equal)


fitHRM.er.fitz = fitHRM(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ER",
                        pi = 'fitzjohn')
plot(fitHRM.er.fitz)

save(fitHRM.er.null, file = './evo.models/nymph-hrm-er-null-antena.RData')
save(fitHRM.er.equal, file = './evo.models/nymph-hrm-er-equal-antena.RData')
save(fitHRM.er.fitz, file = './evo.models/nymph-hrm-er-fitz-antena.RData')

load(file = './evo.models/nymph-hrm-er-null-antena.RData')
load(file = './evo.models/nymph-hrm-er-equal-antena.RData')
load(file = './evo.models/nymph-hrm-er-fitz-antena.RData')


fitHRM.null = fitHRM(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ARD",
                     pi = c(0.5,0,0.5,0,0,0,0,0))
plot(fitHRM.null)


fitHRM.equal = fitHRM(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ARD",
                      pi = 'equal')
plot(fitHRM.equal)


fitHRM.fitz = fitHRM(ladderize(genera_nymph_clock), ant.nymph.sim , model = "ARD",
                     pi = 'fitzjohn')
plot(fitHRM.fitz)

save(fitHRM.null, file = './evo.models/nymph-hrm-null-antena.RData')
save(fitHRM.equal, file = './evo.models/nymph-hrm-equal-antena.RData')
save(fitHRM.fitz, file = './evo.models/nymph-hrm-fitz-antena.RData')

load(file = './evo.models/nymph-hrm-null-antena.RData')
load(file = './evo.models/nymph-hrm-equal-antena.RData')
load(file = './evo.models/nymph-hrm-fitz-antena.RData')


AIC(fit.ER.null, fit.ARD.null,
    fitHRM.er.null, fitHRM.null,
    fit.ER.equal, fit.ARD.equal,
    fitHRM.er.equal, fitHRM.equal,
    fit.ER.fitz, fit.ARD.fitz,
    fitHRM.er.fitz, fitHRM.fitz
)

AIC(fit.ER.equal, fit.ARD.null,
    fitHRM.er.null, fitHRM.null)

AIC(fitHRM.er.null) - AIC(fit.ARD.fitz)

plot(fit.ARD.fitz)
print(fit.ARD.fitz)
anc.ARD.nymph = ancr(fit.ARD.fitz)

pdf(file = './ancestral-tree-nymph-genera.pdf', h = 28, w = 20)
plot(anc.ARD.nymph,  
     args.plotTree = list(lwd = 2, fsize = 1, 
                          ftype = 'i',
                          outline = T,
                          type = 'phylogram',
                          cex = 0.5, offset = 0.6),
     args.nodelabels = list(cex = 0.4),
     args.tiplabels = list(cex = 0.2))
dev.off()

anc.ARD.null = ancr(fit.ARD.null)

pdf(file = './ancestral-tree-null.pdf', h = 28, w = 20)
plot(anc.ARD.null, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()



plot(fit.ARD.corr.null)
anc.ard.corr = ancr(fit.ARD.corr.null)

pdf(file = './ancestral-tree-corr.pdf', h = 28, w = 20)
plot(anc.ard.corr, fsize = 0.7, cex = 0.3,
     ftype = 'i', lwd = 2,
     outline = T)
dev.off()


#### PLOT FULL TREE WITH MISSING DATA ####

name.check(corr.t2,tips)

length(tips$Antena)

antena.full = data.frame(Antena = as.factor(tips$Antena))
rownames(antena.full) = rownames(tips)


name.check(corr.t2,antena.full)

colors<-list(
  setNames(c("black","white", "grey"),c("Expansion","No-data","Straight")))

pdf(file = './full-tree.pdf', h = 10, w = 18)
t = plotFanTree.wTraits(corr.t2, antena.full,
                        colors = colors,
                        spacer = 0,
                        fsize = 0.7, part = 0.5)
dev.off()  

plotTree.datamatrix(corr.t2, antena.full,
                    colors = colors,
                    spacer = 0,
                    fsize = 0.7, part = 0.5)

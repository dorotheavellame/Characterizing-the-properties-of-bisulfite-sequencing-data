#######################################################
### Scripts used to create plots in RRBS data paper ###
### Dorothea Seiler Vellame ###########################
### ds420@exeter.ac.uk ################################
#######################################################


### Figure 1 ##########################################
### RRBS distributions
library(ggplot2)
library(cowplot)
library(reshape2)

load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")

methCols = c(2:126)
covCols = c(127:251)
RRBSMatrixish = allOrdered

## filter out very high coverage, they make the plot too hard to view
RRBSCsub = as.matrix(RRBSMatrixish[,covCols])

x = melt(RRBSCsub)
x = x[complete.cases(x),]
colnames(x) = c("X1", "X2", "value")

plotCovRRBS =
  ggplot(data = x, aes(x = value)) +
  geom_line(aes(group = X2, y = ..count..), stat="density", alpha=0.2) +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  xlab("Read depth") +
  ylab("Number of sites")

RRBSMatrixish[RRBSMatrixish > 200] = NA

RRBSCsub = as.matrix(RRBSMatrixish[,covCols])

x = melt(RRBSCsub)
x = x[complete.cases(x),]
colnames(x) = c("X1", "X2", "value")

plotCovRRBS200 =
  ggplot(data = x, aes(x = value)) +
  geom_line(aes(group = X2, y = ..count..), stat="density", alpha=0.2) +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  xlab("Read depth") +
  ylab("Number of sites")


## do same for the methylation plot
RRBSMsub = as.matrix(RRBSMatrixish[,methCols])/100
x = melt(RRBSMsub[,])

## plotting stats (% in each end 5% of data)
x = x[complete.cases(x),]
nrow(x[x[,3] < 0.0005,])/ nrow(x)
nrow(x[x[,3] > 0.0095,])/ nrow(x)

plotMethRRBS =
  ggplot(data = x, aes(x = value)) +
  geom_line(aes(group = Var2, y = ..count..), stat="density", alpha=0.2) +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  xlab("Proportion of DNA methylation") +
  ylab("Number of sites")

table(plotCovRRBS200$data$value)

save(plotCovRRBS, plotCovRRBS200, plotMethRRBS, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/RRBSDistribution.Rdata")


### Figure 2 ##########################################
### arrayVSRRBSWithFiltering 
load("/mnt/data1/Thea/Simulations/data/ValidationData/longDat.Rdata")

# ## get mean difference between RRBS and array
# mean(longDat$arraylong - longDat$RRBSmethlong)
# sd(longDat$arraylong - longDat$RRBSmethlong)
# 
# ## get mean and sd of sites in each samples covered in the array
# perSample = table(longDat$Sample)
# mean(perSample)
# sd(perSample)

library(ggplot2)
library(cowplot)

rmse = function(m, o){
  sqrt(mean((m - o)^2))
}

methylationCorPlot = function(df,
                              plotCol = 3,
                              colourname = "Coverage",
                              colrange = c(0,125),
                              annotationshift = 0){
  # df must include RRBSmethlong, arraylong, Coverage
  print(ggplot(data = df, aes(x = RRBSmethlong, y = arraylong, col = df[,plotCol])) +
          geom_point() +
          labs(colour = colourname) +
          theme_cowplot(18) +
          xlab("Proportion of methylation \n RRBS") +
          ylab("Proportion of methylation \n Array") +
          scale_color_gradient(low="blue", high="red", limits = colrange, name = "Read depth") +
          labs(fill = "Read depth"))
}

RRBSArrayCorPlot = methylationCorPlot(longDat)


covVals = seq(0, 66, 2)
corFiltP = corFiltRMSE = c()
# corFiltPs = corFiltRMSEs = c()
for (covFilterIt in covVals){
  df = longDat[which(longDat$Coverage > covFilterIt),]
  corFiltP = c(corFiltP, cor(df$arraylong, df$RRBSmethlong))
  
  corFiltRMSE = c(corFiltRMSE, rmse(df$arraylong, df$RRBSmethlong))
  
}

RRBSArraySimilarityDat = data.frame(c(covVals, covVals),
                                    c(corFiltP, corFiltRMSE),
                                    c(rep("Pearson", each = length(corFiltP)),
                                      rep("RMSE", each = length(corFiltP))))
colnames(RRBSArraySimilarityDat) = c("covVals", "corr", "type")

RRBSArrayCorLineDat = RRBSArraySimilarityDat[which(RRBSArraySimilarityDat$type == "Pearson"),]
RRBSArrayRMSEDat = RRBSArraySimilarityDat[which(RRBSArraySimilarityDat$type == "RMSE"),]


RRBSArrayCorLinePlot = ggplot(data = RRBSArrayCorLineDat, aes(x = covVals, y = corr, col = type)) +
  geom_line(cex = 1.5) +
  theme_cowplot(18) +
  xlab("Read depth filter") +
  ylab("Correlation")  +
  ggtitle("") +
  theme(legend.title = element_blank())

RRBSArrayRMSEPlot = ggplot(data = RRBSArrayRMSEDat, aes(x = covVals, y = corr, col = type)) +
  geom_line(cex = 1.5) +
  scale_color_manual(values = c("#1A48C4")) +
  theme_cowplot(18) +
  xlab("Read depth filter") +
  ylab("Error")  +
  ggtitle("") +
  theme(legend.title = element_blank())

save(RRBSArrayCorLinePlot, RRBSArraySimilarityDat, RRBSArrayCorPlot, RRBSArrayRMSEPlot,
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/arrayVSRRBSWithFiltering.Rdata")

### rdCorrDNAmRealSim 
library(ggplot2)
library(cowplot)
library(reshape2)

load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")

DNAm = melt(allOrdered[,2:126])
cov = melt(allOrdered[,127:251])
flatsamples = cbind.data.frame(DNAm, cov[,2])
flatsamples = flatsamples[complete.cases(flatsamples),]
flatsamples = flatsamples[which(flatsamples[,3] > 10),]
flatsamples$decile = cut(flatsamples$value, seq(0,100,10), include.lowest = T)

save(flatsamples, file = "/mnt/data1/Thea/Simulations/data/flatsamplesForrdCorrDNAmRealSim.Rdata")


load("/mnt/data1/Thea/Simulations/data/flatsamplesForrdCorrDNAmRealSim.Rdata")
nPerSample = 100
subflat = data.frame(matrix(nrow = nPerSample*10, ncol = 4))
declev = levels(flatsamples$decile)
for (i in 1:10){
  temp = flatsamples[which(flatsamples$decile == declev[i]),]
  subflat[((i-1)*nPerSample+1):(i*nPerSample),] = temp[sample(rownames(temp), nPerSample),]
}
save(flatsamples, subflat, file = "/mnt/data1/Thea/Simulations/data/flatsamplesForrdCorrDNAmRealSim.Rdata")

load("/mnt/data1/Thea/Simulations/data/flatsamplesForrdCorrDNAmRealSim.Rdata")
colnames(subflat) = c("sample", "DNAm", "cov", "decile")

simulateData = function(rd,
                        muVec,
                        r = 1.5){
  obs = sapply(muVec, function(mu){return(rbinom(1, rd, mu)/rd)})
  return(obs)
}

simDat = sapply(c(1:50), simulateData, subflat$DNAm/100)

colnames(simDat) = 1:50
rmse = function(m, o){
  sqrt(mean((m - o)^2))
}
correlations = spcor = RMSe =c()
for (i in 1:50){
  correlations[i] = cor(subflat$DNAm/100, simDat[,i], method = "pearson")
  RMSe[i] = rmse(subflat$DNAm/100, simDat[,i])
}

plotDatCor = cbind.data.frame(cor =  correlations, rd = 1:50, error = RMSe)
rdSimCorPlot =
  ggplot(plotDatCor, aes(x = rd, y = cor)) +
  geom_line(size = 2, col = "#F8766D") +
  theme_cowplot(18) +
  labs(x = "Read depth", y = "Correlation")

rdSimErrorPlot =
  ggplot(plotDatCor, aes(x = rd, y = error)) +
  geom_line(size = 2, col = "#1A48C4") +
  theme_cowplot(18) +
  labs(x = "Read depth", y = "Error")

save(rdSimCorPlot, rdSimErrorPlot, file = "/mnt/data1/Thea/Simulations/data/rdCorrDNAmRealSimPlot.Rdata")


### simulatedPrecision 
simulateData = function(rd,
                        r = 1.5,
                        mu = 0.05,
                        nPerm = 1000){
  for(i in 1:nPerm){
    obs1<-rbinom(nPerm, rd, mu)/rd
  }
  return(cbind.data.frame(DNAmPred = obs1, rd = rd, DNAmmu = mu))
}

x = lapply(c(1,2,3,4,5,10,15,20,25,30,35,40,45,50),simulateData)

y = x[[1]]
for(i in 2:length(x)){
  y = rbind.data.frame(y,x[[i]])
}

y$rd = as.factor(y$rd)

library(ggplot2)
library(cowplot)
simulatedPrecision = ggplot(y, aes(x = rd, y = DNAmPred, fill = rd)) +
  geom_violin() +
  theme_cowplot(18) +
  geom_hline(yintercept = 0.5, col = "red", linetype = "dashed") +
  theme(legend.position = "none") +
  labs(y = "DNAm", x = "Read depth")


for(i in 1:length(seq(0,1,0.1))){
  x[[i]] = lapply(1:10,simulateData, mu = seq(0,1,0.1)[i], nPerm = 1000)
}

z = data.frame(matrix(ncol = 3, nrow = 0))
for(i in 1:11){
  for(j in 1:10){
    z = rbind.data.frame(z,x[[i]][[j]])
  }
}

z$rd = as.factor(z$rd)
countDataPlot = ggplot(z, aes(x = rd, y = DNAmPred, fill = rd)) +
  geom_violin() +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  labs(y = "DNAm", x = "Read depth")
countDataPlot

save(simulatedPrecision, countDataPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/simulatedPrecision.Rdata")

### proportion in extreme DNAm 
# load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
# 
# RRBSmat = allOrdered
# 
# inExtreme = function(rdnvec, RRBSmat){
#   methCol = RRBSmat[,rdnvec[2] +1]
#   covCol = RRBSmat[,rdnvec[2] +126]
#   
#   methRd = methCol[which(covCol == rdnvec[1])]
#   
#   lex = sum(methRd<5)
#   uex = sum(methRd>95)
#   
#   return(c(lex+uex, length(methRd)))
# }
# 
# rdn = expand.grid(c(1,seq(5,50,5)), 1:125)
# propInEx = t(apply(rdn, 1, inExtreme, RRBSmat))
# 
# Dat = data.frame(extreme = propInEx[,1], total = propInEx[,2], rd = rdn[,1], nSample = rdn[,2])
# 
# save(Dat, file = "/mnt/data1/Thea/Simulations/data/propInExtremes.Rdata")

load("/mnt/data1/Thea/Simulations/data/propInExtremes.Rdata")

Dat$percentage = Dat$extreme/Dat$total*100

# ## stats
# mean(Dat$percentage[which(Dat$rd == 5)])
# mean(Dat$percentage[which(Dat$rd == 50)])
# sd(Dat$percentage[which(Dat$rd == 5)])
# sd(Dat$percentage[which(Dat$rd == 50)])


library(ggplot2)
library(cowplot)

extremesPlot = ggplot(Dat, aes(x = as.factor(rd), y = percentage, fill = as.factor(rd))) +
  geom_boxplot() +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  labs(x = "Read depth", y = "DNAm points in extreme (%)")

save(extremesPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/extremesPlot.Rdata")

### Figure 3 ##########################################
### overlapBetweenSampleCpGs 
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
library(ggplot2)
library(cowplot)

cov = allOrdered[,127:ncol(allOrdered)]
xMax = 25
nPoss = seq(5,125,5)

plotDat = as.data.frame(matrix(nrow = xMax*length(nPoss), ncol = 2))

i = 1
for (x in 1:xMax){ # number of times to repeate sampling
  for (n in nPoss){ # number of samples used in intersection calculation
    nSelected = sample(1:ncol(cov), n)
    plotDat[i,1] = n
    plotDat[i,2] = sum(complete.cases(cov[,nSelected]))/ nrow(cov[which(rowSums(is.na(cov[,nSelected])) != n),])
    i = i + 1
  }
}

colnames(plotDat) = c("n","prop")
plotDat$n = as.factor(plotDat$n)

overlapBetweenSamplesPlotDat = plotDat
overlapBetweenSamplesPlot = ggplot(aes(x = n, y = prop), data = plotDat) +
  geom_boxplot() +
  theme_cowplot(18) +
  xlab("Number of samples compared") +
  ylab("Proportion of DNAm sites present \nin all samples compared") +
  ylim(c(0.2,0.9))

## save the plot data incase plot needs adjusting
save(overlapBetweenSamplesPlot, overlapBetweenSamplesPlotDat,
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/OverlapBetweenSampleCpGs.Rdata")

### corrHistBetweenSamples 
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
library(ggplot2)
library(cowplot)
library(reshape)

datasub = allOrdered[,127:251]
for (i in 1:ncol(datasub)){
  datasub[is.na(datasub[,i]),i] = 0
}

cordat = cor(datasub, method = "spearman")
## remove upper traingle and diag and unlist
for(i in 1:nrow(cordat)){
  for (j in 1:ncol(cordat)){
    if (j >= i){
      cordat[i,j] = NA
    }
  }
}

cordat = melt(cordat)
cordat = cordat[complete.cases(cordat[,3]),]

rdCorHist = ggplot(data = cordat, aes(x = value)) +
  geom_histogram(bins = 50) +
  theme_cowplot(18) +
  xlab("Correlation") +
  ylab("Number of pairs") +
  xlim(c(0,1))

save(rdCorHist,
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/corrHistBetweenSamples.Rdata")

### dataRemaining 
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
sub = allOrdered[,c(127:251)]

dataRemaining = function(dataCol, rdFilt){
  dataCol[which(dataCol<rdFilt)] = NA
  return(sum(!is.na(dataCol))/length(dataCol))
}

rdRange = 1:150

propRemaining = data.frame(matrix(nrow = length(rdRange), ncol = 125))
for (rd in rdRange){
  propRemaining[rd,] = apply(sub,2, dataRemaining, rd)
}
colnames(propRemaining) = colnames(sub)

plotDat = data.frame(value = apply(propRemaining, 1, mean),
                     rd = rdRange)

propRemainPlot = ggplot(plotDat, aes(y = value, x = rd)) +
  geom_line(size = 1.5) +
  theme_cowplot(18) +
  labs(x = "Read depth filter", y = "Proportion of data remaining") +
  ylim(c(0,1))


save(propRemainPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/dataRemaining.Rdata")

### rdCorCol 
library(ggplot2)
library(cowplot)
library(reshape)
library(viridis)

## load data
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
sub = allOrdered[sample(1:nrow(allOrdered),1000),c(127:251)]

nSamplesPerSite = apply(sub, 1, function(x){return(125-sum(is.na(x)))})

covDat = melt(sub)
covDat$nSamplesPerSite = nSamplesPerSite


temp = covDat[which(covDat$variable == "A17_cov"),]

temp$A18 = covDat[which(covDat$variable == "A18_cov"),2]

rdCorColPlot = ggplot(temp, aes(x = value, y = A18, col = nSamplesPerSite)) +
  geom_point(alpha = 0.8) +
  theme_cowplot(18) +
  xlim(c(0,200)) +
  ylim(c(0,200)) +
  scale_color_viridis(option = "D", direction = -1) +
  labs(col = "Samples\nper site", x = "Read depth in sample 1", y = "Read depth in sample 2")

save(rdCorColPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/rdCorCol.Rdata")

### annotationPlot 
## load packages
library(annotatr)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(viridis)

## load data
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")

## create GRanges for all RRBS data
## separate chr from range
location = unlist(strsplit(allOrdered[,1], "_"))
chr = location[seq(1,length(location), 2)]
range = as.numeric(location[seq(2,length(location), 2)])

grALL = GRanges(
  seqnames = Rle(chr),
  ranges = IRanges(range, end = range+1)
)

## create GRanges for RRBS where sites are in all samples
sub = allOrdered[complete.cases(allOrdered),]
location = unlist(strsplit(sub[,1], "_"))
chr = location[seq(1,length(location), 2)]
range = as.numeric(location[seq(2,length(location), 2)])

gr125 = GRanges(
  seqnames = Rle(chr),
  ranges = IRanges(range, end = range+1)
)

## make annotations for grALL and gr125
annotcpg = build_annotations(genome = "mm10", annotations = builtin_annotations()[grep("mm10_cpg", builtin_annotations())])
annotgene = build_annotations(genome = "mm10", annotations = builtin_annotations()[grep("mm10_genes", builtin_annotations())])

aALL = list(data.frame(annotate_regions(grALL, annotations = annotcpg  , minoverlap = 1L, ignore.strand = TRUE, quiet = FALSE)),
            data.frame(annotate_regions(grALL, annotations = annotgene  , minoverlap = 1L, ignore.strand = TRUE, quiet = FALSE)))

a125 = list(data.frame(annotate_regions(gr125, annotations = annotcpg  , minoverlap = 1L, ignore.strand = TRUE, quiet = FALSE)),
            data.frame(annotate_regions(gr125, annotations = annotgene  , minoverlap = 1L, ignore.strand = TRUE, quiet = FALSE)))

pieDatcpg = rbind.data.frame(cbind.data.frame(table(aALL[[1]]$annot.type), group = "ALLcpg"),
                             cbind.data.frame(table(a125[[1]]$annot.type), group = "125cpg"))

pieDatgene = rbind.data.frame(cbind.data.frame(table(aALL[[2]]$annot.type), group = "ALLgene"),
                              cbind.data.frame(table(a125[[2]]$annot.type), group = "125gene"))

annotationPlot = ggplot(pieDatcpg, aes(x = group, y = Freq, fill = Var1)) +
  geom_bar(position="fill", stat="identity") +
  theme_cowplot(18) +
  scale_fill_viridis(discrete = T, name = "CpG location", labels = c("Inter", "Islands", "Shelves", "Shores")) +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank()) +
  scale_x_discrete(labels = c("All sites", "Sites in all \nsamples"))

save(pieDatcpg, annotationPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/annotationPlot.Rdata")
### Figure 4 ##########################################
library(ggplot2)
library(cowplot)
library(colortools)

simulateData = function(n1 = 20,
                        n2 = 20,
                        r = 1.5,
                        mu = 25,
                        meanDiff = 20,
                        rd.mu = 20,
                        nPerm = 10000){
  sim = matrix(data = NA, ncol = nPerm, nrow = 1)
  mu1<-mu + meanDiff
  mu2<-mu
  for(j in 1:nPerm){
    ## sample read depth from nbinomial distribution with given mean
    rd1<-rnbinom(n1,r,mu=rd.mu)
    rd2<-rnbinom(n2,r,mu=rd.mu)
    while (sum(rd1==0)>0){
      k<-which(rd1==0)
      rd1[k]<-rnbinom(length(k),r,mu=rd.mu)
    }
    while (sum(rd2==0)>0){
      k<-which(rd2==0)
      rd2[k]<-rnbinom(length(k),r,mu=rd.mu)
    }
    ## use binomial distibution to sample from read depth number of methylated reads and calc DNAm value
    obs1<-rep(NA, n1)
    for(i in 1:n1){
      obs1[i]<-rbinom(1, rd1[i], mu1/100)/rd1[i]
    }
    obs2<-rep(NA, n2)
    for(i in 1:n2){
      obs2[i]<-rbinom(1, rd2[i], mu2/100)/rd2[i]
    }
    sim[j]<-t.test(obs1, obs2)$p.value
  }
  return(sim)
}

nPermUse = 10000

## simulate read depth changes
rdPoss = seq(10,100,1)
rdS1 = rdS2 = rdS3 = matrix(ncol = length(rdPoss), nrow = nPermUse)

for(rdIndex in 1:length(rdPoss)){
  rdS1[,rdIndex] = simulateData(rd.mu = rdPoss[rdIndex], nPerm = nPermUse,
                                n1 = 30,
                                n2 = 30,
                                r = 1.5,
                                mu = 25,
                                meanDiff = 20)
  rdS2[,rdIndex] = simulateData(rd.mu = rdPoss[rdIndex], nPerm = nPermUse,
                                n1 = 30,
                                n2 = 30,
                                r = 1.5,
                                mu = 25,
                                meanDiff = 10)
  rdS3[,rdIndex] = simulateData(rd.mu = rdPoss[rdIndex], nPerm = nPermUse,
                                n1 = 30,
                                n2 = 30,
                                r = 1.5,
                                mu = 25,
                                meanDiff = 5)
}

## reformat for plotting
pS1 = colSums(rdS1[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS2 = colSums(rdS2[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS3 = colSums(rdS3[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100

rdPlotDat = rbind.data.frame(cbind.data.frame(dat = pS1, group = "i", rd = rdPoss),
                             cbind.data.frame(dat = pS2, group = "ii", rd = rdPoss),
                             cbind.data.frame(dat = pS3, group = "iii", rd = rdPoss))


rdPlot = ggplot(rdPlotDat, aes(x = rd, y = dat, col = group)) +
  geom_line(size = 1.5) +
  ylim(0,100) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  xlab("Read depth") +
  ylab("Power (%)") +
  theme_cowplot(18) +
  labs(col = element_blank()) +
  scale_color_manual(values = analogous("#1AC496"))



### sample size and sample size balance and power ##
# set default parameters
# n1 = n2 = 20
r = 1.5
mu = 25 # mean methylation of sites
meanDiff = 20
rd.mu = 20 # mean read depth
nPerm<-nPermUse

meanNPoss = 5:100 # possible sample sizes
n1BalancePoss = c(0.5, 0.6, 0.7, 0.8)
n2BalancePoss = c(0.5, 0.4, 0.3, 0.2)

# enough rows for all balances
sim <- matrix(data = NA, ncol = nPerm , nrow = length(meanNPoss)*length(n1BalancePoss))

for (nBalanceIndex in 1:length(n1BalancePoss)){
  
  mu1<-mu + meanDiff
  mu2<-mu
  
  for (nIndex in 1:length(meanNPoss)){
    
    n1 = round(meanNPoss[nIndex]*n1BalancePoss[nBalanceIndex]*2)
    n2 = round(meanNPoss[nIndex]*n2BalancePoss[nBalanceIndex]*2)
    
    for(j in 1:nPerm){
      ## sample read depth from nbinomial distribution with given mean
      rd1<-rnbinom(n1,r,mu=rd.mu)
      rd2<-rnbinom(n2,r,mu=rd.mu)
      while (sum(rd1==0)>0){
        k<-which(rd1==0)
        rd1[k]<-rnbinom(length(k),r,mu=rd.mu)
      }
      while (sum(rd2==0)>0){
        k<-which(rd2==0)
        rd2[k]<-rnbinom(length(k),r,mu=rd.mu)
      }
      ## use binomial distibution to sample from read depth number of methylated reads and calc DNAm value
      obs1<-rep(NA, n1)
      for(i in 1:n1){
        obs1[i]<-rbinom(1, rd1[i], mu1/100)/rd1[i]
      }
      obs2<-rep(NA, n2)
      for(i in 1:n2){
        obs2[i]<-rbinom(1, rd2[i], mu2/100)/rd2[i]
      }
      
      sim[(nBalanceIndex-1)*length(meanNPoss)+nIndex,j]<-t.test(obs1, obs2)$p.value
      
    }
  }
}

balanceStatus = c(rep("50/50",length(meanNPoss)),
                  rep("60/40",length(meanNPoss)),
                  rep("70/30",length(meanNPoss)),
                  rep("80/20",length(meanNPoss)))

# reformat data for plotting
sampleSizeDat = data.frame(c(meanNPoss,meanNPoss,meanNPoss,meanNPoss)*2,
                           rowSums(sim[,1:nPerm] < (0.05/nPerm))/nPerm*100,
                           balanceStatus)
colnames(sampleSizeDat) = c("meanNPoss", "power", "balanceGroup")

sampleSizePlot =
  ggplot(sampleSizeDat, aes(x = meanNPoss, y = power, col = balanceGroup)) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  ylim(0,100) +
  xlab("Total sample size") +
  ylab("Power (%)") +
  labs(col = expression(atop("Sample", paste("split")))) +
  theme_cowplot(18) +
  scale_color_manual(values = unique(c(analogous("#961AC4"),analogous("#C41A9D"))))


## simulate meanDiff changes
mdPoss = seq(1,50,0.5)
mdS1 = mdS2 = mdS3 = mdS4 = matrix(ncol = length(mdPoss), nrow = nPermUse)
mdS5 = mdS6 = matrix(ncol = length(mdPoss), nrow = nPermUse)
for(mdIndex in 1:length(mdPoss)){
mdS1[,mdIndex] = simulateData(rd.mu = 25, nPerm = nPermUse,
                              n1 = 20,
                              n2 = 20,
                              r = 1.5,
                              mu = 25,
                              meanDiff = mdPoss[mdIndex])
mdS2[,mdIndex] = simulateData(rd.mu = 25, nPerm = nPermUse,
                              n1 = 10,
                              n2 = 10,
                              r = 1.5,
                              mu = 25,
                              meanDiff = mdPoss[mdIndex])
  mdS3[,mdIndex] = simulateData(rd.mu = 25, nPerm = nPermUse,
                                n1 = 5,
                                n2 = 5,
                                r = 1.5,
                                mu = 25,
                                meanDiff = mdPoss[mdIndex])
## final added after reviewers comments
  mdS4[,mdIndex] = simulateData(rd.mu = 25, nPerm = nPermUse,
                                n1 = 50,
                                n2 = 50,
                                r = 1.5,
                                mu = 25,
                                meanDiff = mdPoss[mdIndex])
  
  mdS5[,mdIndex] = simulateData(rd.mu = 25, nPerm = nPermUse,
                                n1 = 100,
                                n2 = 100,
                                r = 1.5,
                                mu = 25,
                                meanDiff = mdPoss[mdIndex])
  
  mdS6[,mdIndex] = simulateData(rd.mu = 25, nPerm = nPermUse,
                                n1 = 500,
                                n2 = 500,
                                r = 1.5,
                                mu = 25,
                                meanDiff = mdPoss[mdIndex])
  
  
}

## reformat for plotting
pS1 = colSums(mdS1[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS2 = colSums(mdS2[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS3 = colSums(mdS3[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS4 = colSums(mdS4[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS5 = colSums(mdS5[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS6 = colSums(mdS6[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100

mdPlotDat = rbind.data.frame(cbind.data.frame(dat = pS3, group = "i", md = mdPoss/100),
                             cbind.data.frame(dat = pS2, group = "ii", md = mdPoss/100),
                             cbind.data.frame(dat = pS1, group = "iii", md = mdPoss/100),
                             cbind.data.frame(dat = pS4, group = "iv", md = mdPoss/100),
                             cbind.data.frame(dat = pS5, group = "v", md = mdPoss/100),
                             cbind.data.frame(dat = pS6, group = "vi", md = mdPoss/100))

mdPlot = ggplot(mdPlotDat, aes(x = md, y = dat, col = group)) +
  geom_line(size = 1.5) +
  ylim(0,100) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  xlab("Mean DNAm difference") +
  ylab("Power (%)") +
  theme_cowplot(18) +
  labs(col = element_blank())

simulateData = function(n1 = 20,
                        n2 = 20,
                        r = 1.5,
                        mu = 25,
                        meanDiff = 20,
                        rd.mu = 20,
                        nPerm = 10000,
                        varCalc = F){
  sim = varSim = matrix(data = NA, ncol = nPerm, nrow = 1)
  mu1<-mu + meanDiff/2
  mu2<-mu - meanDiff/2
  for(j in 1:nPerm){
    ## sample read depth from nbinomial distribution with given mean
    rd1<-rnbinom(n1,r,mu=rd.mu)
    rd2<-rnbinom(n2,r,mu=rd.mu)
    while (sum(rd1==0)>0){
      k<-which(rd1==0)
      rd1[k]<-rnbinom(length(k),r,mu=rd.mu)
    }
    while (sum(rd2==0)>0){
      k<-which(rd2==0)
      rd2[k]<-rnbinom(length(k),r,mu=rd.mu)
    }
    ## use binomial distibution to sample from read depth number of methylated reads and calc DNAm value
    obs1<-rep(NA, n1)
    for(i in 1:n1){
      obs1[i]<-rbinom(1, rd1[i], mu1/100)/rd1[i]
    }
    obs2<-rep(NA, n2)
    for(i in 1:n2){
      obs2[i]<-rbinom(1, rd2[i], mu2/100)/rd2[i]
    }
    if (varCalc == T){
      varSim[j] = var(obs1)
    }
    sim[j]<-t.test(obs1, obs2)$p.value
  }
  if(varCalc == T){
    return(varSim)
  }
  return(sim)
}


## simulate meanDiff changes
muPoss = seq(5, 95, 1)
muS1 = muS2 = muS3 =
  vS1 = vS2 = vS3 = matrix(ncol = length(muPoss), nrow = nPermUse)

muS1 = sapply(muPoss, simulateData,
              rd.mu = 50,
              nPerm = nPermUse,
              n1 = 80,
              n2 = 80,
              r = 1.5,
              meanDiff = 5)
muS2 = sapply(muPoss, simulateData,
              rd.mu = 30,
              nPerm = nPermUse,
              n1 = 80,
              n2 = 80,
              r = 1.5,
              meanDiff = 5)
muS3 = sapply(muPoss, simulateData,
              rd.mu = 10,
              nPerm = nPermUse,
              n1 = 80,
              n2 = 80,
              r = 1.5,
              meanDiff = 5)
vS1 = sapply(muPoss, simulateData,rd.mu = 50, nPerm = nPermUse,
             n1 = 80,
             n2 = 80,
             r = 1.5,
             meanDiff = 5,
             varCalc = T)
vS2 = sapply(muPoss, simulateData, rd.mu = 30, nPerm = nPermUse,
             n1 = 80,
             n2 = 80,
             r = 1.5,
             meanDiff = 5,
             varCalc = T)
vS3 = sapply(muPoss, simulateData,rd.mu = 10, nPerm = nPermUse,
             n1 = 80,
             n2 = 80,
             r = 1.5,
             meanDiff = 5,
             varCalc = T)

## reformat for plotting
pS1 = colSums(muS1[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS2 = colSums(muS2[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100
pS3 = colSums(muS3[1:nPermUse,] < (0.05/nPermUse))/nPermUse*100

v1 = apply(vS1, 2, mean)
v2 = apply(vS2, 2, mean)
v3 = apply(vS3, 2, mean)

muPlotDat = rbind.data.frame(cbind.data.frame(dat = pS1, var = v1, group = "i", mu = muPoss),
                             cbind.data.frame(dat = pS2, var = v2, group = "ii", mu = muPoss),
                             cbind.data.frame(dat = pS3, var = v3, group = "iii", mu = muPoss))

muPlot = ggplot(muPlotDat, aes(x = mu/100, y = dat, col = group)) +
  geom_line(size = 1.5) +
  ylim(0,100) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  theme_grey(base_size = 18) +
  xlab("Mean DNAm") +
  ylab("Power (%)") +
  theme_cowplot(18) +
  labs(col = element_blank()) +
  scale_color_manual(values = analogous("#9DC41A"))

varPlot = ggplot(muPlotDat, aes(x = mu/100, y = var, col = group)) +
  geom_line(size = 1.5) +
  theme_grey(base_size = 18) +
  xlab("Mean DNAm") +
  ylab("Variance") +
  theme_cowplot(18) +
  labs(col = element_blank()) +
  scale_color_manual(values = analogous("#9DC41A"))

save(mdPlot, muPlot, rdPlot, sampleSizePlot, varPlot,
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/simPlotData.Rdata")



### Figure 5 ##########################################
### cumulitiveOverRD
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
RRBSmat = allOrdered


inExtreme = function(rdnvec, RRBSmat){
  covCol = RRBSmat[,rdnvec[2] +126]
  covRd = covCol[which(covCol >= rdnvec[1])]
  return(c(length(covRd)))
}

rdn = expand.grid(c(1,seq(5,300,5)), 1:125)
covOverRd = t(apply(rdn, 1, inExtreme, RRBSmat))

totalPoints = t(apply(cbind(rep(0, 125), 1:125), 1, inExtreme, RRBSmat))
tpSum = sum(totalPoints)

Dat = data.frame(extreme = t(covOverRd), rd = rdn[,1], nSample = rdn[,2])

plotDat = data.frame(matrix(nrow=length(c(1,seq(5,300,5))), ncol = 2))
colnames(plotDat) = c("totOver", "rd")
plotDat$rd = c(1,seq(5,300,5))
for (i in 1:length(plotDat$rd)){
  plotDat$totOver[i] = sum(Dat[which(Dat$rd == plotDat$rd[i]),"extreme"])
}

library(cowplot)
library(ggplot2)

plotDat$tot = plotDat$totOver/tpSum*100

minRDPlot =
ggplot(plotDat, aes(x = rd, y = tot))+
  geom_line() +
  theme_cowplot(18) +
  labs(x = "Minimum read depth", y = "Data remaining (%)")
save(minRDPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/minRDPlot.Rdata")

### simQQ 
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
RRBSMatrix = allOrdered

source("/mnt/data1/Thea/Simulations/RScripts/POWEREDBiSeqFunction.R")

## set filtering
rdFilt = 1
meanRD = mean(as.matrix(RRBSMatrix[,127:251]), na.rm = T)

## filter data
rdFiltMatrix = function(rowInd, rdFilt, covCols, methCols, nSamplesNeededBoth, RRBSIn){
  if(sum(!is.na(RRBSIn[rowInd, covCols]))< nSamplesNeededBoth){
    RRBSIn[rowInd, 1] = NA
    return(RRBSIn[rowInd,])
  }
  RRBSIn[rowInd, na.omit(methCols[RRBSIn[rowInd, covCols] < rdFilt])] = NA
  RRBSIn[rowInd, na.omit(covCols[RRBSIn[rowInd, covCols] < rdFilt])] = NA
  if(sum(!is.na(RRBSIn[rowInd, covCols]))< nSamplesNeededBoth){
    RRBSIn[rowInd, 1] = NA
    return(RRBSIn[rowInd,])
  }
  return(RRBSIn[rowInd,])
}

rows = 1:10000
z = t(sapply(rows ,FUN = rdFiltMatrix, rdFilt = rdFilt, covCols = 3, methCols = 2, nSamplesNeededBoth = 1, RRBSIn = RRBSMatrix[rows,c(1,2,127)]))
RRBSFilt = z[!is.na(z[,1]),]

priorCalc = function(RRBSm){
  if(is.null(ncol(RRBSm))){
    RRBSm = cbind(RRBSm, RRBSm)
  }
  
  sumCalc = function(RRBSmRow){
    g1 = sum(RRBSmRow < 5, na.rm = T)
    g3 = sum(RRBSmRow > 95, na.rm = T)
    tot = sum(!is.na(RRBSmRow))
    return(c(g1, g3, tot))
  }
  
  propDat = colSums(t(apply(RRBSm, 1, sumCalc)))
  priors = c(propDat[1]/propDat[3],
             (propDat[3]-(propDat[1]+propDat[2]))/propDat[3],
             propDat[2]/propDat[3])
  
  return(priors)
}


muDistrib = list(function(){runif(1, min = 0, max = 5)},
                 function(){runif(1, min = 5, max = 95)},
                 function(){runif(1, min = 95, max = 100)})

simFuncQQ = function(nSamples, rdFilt, rdSim, datForR, nSamplesNeededBoth, prior){
  r = parameterChecksFunction(matrix(data = datForR[1:round(length(datForR), digits = -1)],
                                     ncol = 10, 
                                     nrow = round(length(datForR), digits = -1)/10))[4]
  
  rd1 = rnbinom(nSamplesNeededBoth, r, mu = rdSim)
  
  while(any(rd1<rdFilt)){
    rd1[rd1<rdFilt] = rnbinom(sum(rd1<rdFilt), r, mu = rdSim) 
  }
  
  if(nSamplesNeededBoth<nSamples){
    rd1new = rnbinom(nSamples - nSamplesNeededBoth, r, mu = rdSim)
    rd1new = rd1new[rd1new>rdFilt]
    rd1 = c(rd1, rd1new)
  }
  
  ## use binomial distibution to sample from read depth number of methylated reads and calc DNAm value
  obs1<-rep(NA, length(rd1))
  for(i in 1:length(rd1)){
    obs1[i]<-rbinom(1, rd1[i], muDistrib[[sample(3, size = 1, prob = prior)]]()/100)/rd1[i]
  }
  return(cbind(obs1, rd1))
}

prior = priorCalc(RRBSFilt[,2])

simDat = simFuncQQ(nSamples = nrow(RRBSFilt), rdFilt = rdFilt, rdSim = meanRD, nSamplesNeededBoth = nrow(RRBSFilt), prior = prior, 
                   datForR = unlist(RRBSFilt[,3]))

plotDat = data.frame(rdSim = sort(simDat[,2]), DNAmSim = sort(simDat[,1]),
                     rdReal = sort(unlist(RRBSFilt[,3])), DNAmReal = sort(unlist(RRBSFilt[,2]))/100)

library(ggplot2)
library(cowplot)

DNAmQQ =
  ggplot(data = plotDat, aes(x = DNAmReal, y = DNAmSim)) +
  geom_line(size = 1.5) +
  theme_cowplot(18) +
  geom_abline() +
  labs(x = "True DNAm", y = "Simulated DNAm")

RDQQFull =
  ggplot(data = plotDat, aes(x = rdReal, y = rdSim)) +
  geom_point() +
  theme_cowplot(18) +
  geom_abline() +
  labs(x = "True RD", y = "Simulated RD")

save(DNAmQQ, RDQQFull, 
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/simQQ.Rdata")


### plot diff in power at rd 1-75 at various nsamplesneeded
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
source("/mnt/data1/Thea/Simulations/RScripts/POWEREDBiSeqFunction.R")

powernIN = function(x){t(sapply(c(1,75), POWEREDBiSeq,
                                allOrdered, 
                                meanDiff = 0.06, 
                                pheno = FALSE, 
                                nSampleNeeded = x, 
                                optimalSearchNPerm = 40000))}

outList = lapply(seq(2,60,2), powernIN)

save(outList, file ="/mnt/data1/Thea/Simulations/data/fromISCA/diffComparewithNSN1.RData")

library(ggplot2)
library(cowplot)
library(reshape2)

load("/mnt/data1/Thea/Simulations/data/fromISCA/diffComparewithNSN1.RData")
pDiff1 = c()
for (i in 1:length(outList)){
  pDiff1 = c(pDiff1, unlist(outList[[i]][2,2]) -  unlist(outList[[i]][1,2]))
}

pDiff = cbind.data.frame(pDiff1, nSamples = seq(2,60,2))

powerDiffPlot = ggplot(pDiff, aes(x = nSamples, y = pDiff1)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "Number of samples needed", y = "Difference in power")

save(powerDiffPlot, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/diffInPower1_75rd.Rdata")


### how tool power varies
library(ggplot2)
library(cowplot)

## done in ISCA, not knight ##
## use .sh script:
# SBATCH --array=0-9
# nPerm=(2 30 60)
# Rscript subScript.R ${nPerm[${SLURM_ARRAY_TASK_ID}]}
# 
# subScript.R :
# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# load("/gpfs/ts0/projects/Research_Project-191406/Thea/DebugRRBS/data/RRBSmergedThatCorrelate.Rdata")
# source("/gpfs/ts0/projects/Research_Project-191406/Thea/DebugRRBS/RScripts/POWEREDBiSeqFunction.R")
# x = replicate(40,t(POWEREDBiSeq(15,
#                                 allOrdered,
#                                 meanDiff = 0.06,
#                                 pheno = FALSE,
#                                 nSampleNeeded = 60,
#                                 optimalSearchNPerm = 40000)))
# save(x, file = paste("/gpfs/ts0/projects/Research_Project-191406/Thea/DebugRRBS/data/powerOutHistDat_",args[1],".Rdata", sep = ""))

n = c(2,30,60)
powerOut = data.frame(matrix(nrow = 0, ncol = 8))
for (i in 1:3){
  load(paste("/mnt/data1/Thea/Simulations/data/fromISCA/powerOutnSamplesNeeded_", n[i], ".Rdata", sep = ""))
  x = data.frame(x)
  powerOut = rbind.data.frame(powerOut2, data.frame(rd = unlist(x[,1]),
                                                    power = unlist(x[,2]),
                                                    minPower = unlist(x[,3]),  
                                                    maxPower = unlist(x[,4]),     
                                                    bonferroniP = unlist(x[,5]),
                                                    nCpGTestedProp = unlist(x[,6]),
                                                    nCpGTested = unlist(x[,7]),
                                                    nSampleNeeded = n[i]))
}

powerWithPOWEREDBiSeq =
  ggplot(powerOut, aes(x = rd, y = power, col = as.factor(nSampleNeeded), shape = as.factor(nSampleNeeded))) +
  geom_line(size = 1) +
  theme_cowplot(18) +
  geom_point(size = 2) +
  labs(x = "Read depth", y = "Power (%)", col = "Samples\nNeeded", shape = "Samples\nNeeded")

save(powerWithPOWEREDBiSeq, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/powerWithPOWEREDBiSeq.Rdata")


























### Supplementary Figure 1 and 2 ######################
## plot basic distributions of the subsets of data
## read in the data
load("/mnt/data1/Thea/Simulations/data/ValidationData/arrayRRBSComparisonData.Rdata")
library(ggplot2)
library(cowplot)
library(reshape)

## make data plotable
arrayplot = as.data.frame(as.numeric(unlist(c(arraysub))))
methplot = as.data.frame(as.numeric(unlist(c(RRBSmeth))))
colnames(arrayplot) = colnames(methplot) = "data"

arrayplot = melt(arraysub)
methplot = melt(RRBSmeth)
colnames(methplot) = c("X1", "X2", "value")

arrayRRBSDistribDat = data.frame(sample = c(paste(arrayplot$variable, "array", sep = ""), methplot$X2),
                                 meth = c(arrayplot$value, methplot$value),
                                 type = rep(c("array", "RRBS"), each = nrow(methplot)))

arrayRRBSDistribPlot = ggplot(data = arrayRRBSDistribDat, aes(x = meth, col = type)) +
  geom_line(aes(group = sample, y = ..count..), stat="density", alpha = 0.3) +
  theme_cowplot(18) +
  theme(legend.title = element_blank()) +
  xlab("Proportion of DNAm") +
  ylab("Number of sites") +
  scale_color_manual(values=c('#2a18a2','#90a218')) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_y_continuous(labels = scales::label_number_auto())


load("/mnt/data1/Thea/Simulations/data/ValidationData/longDat.Rdata")
longDat$residuals = longDat$RRBSmethlong - longDat$arraylong

ResidualsRRBSArrayPlot = ggplot(data = longDat, aes(x = Coverage, y = residuals)) +
  geom_point(alpha = 0.3) +
  theme_cowplot(18) +
  geom_hline(yintercept = 0, color = "red" ,linetype = "dashed") +
  xlab("Read depth") +
  ylab("RRBS DNAm - array DNAm")

ggplot(data = longDat, aes(x = as.factor(Coverage), y = residuals)) +
  geom_boxplot()+
  theme_cowplot(18) +
  geom_hline(yintercept = 0, color = "red" ,linetype = "dashed") +
  xlab("Read depth") +
  ylab("RRBS DNAm - array DNAm") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.x.bottom = element_text(size=10))

## use the cut function to bin the data
longDat$cutDat = cut(longDat$Coverage, c(seq(0,15,3),seq(25, 45, 10), 60, 90, 124))

ResidualsRRBSArrayPlot = ggplot(data = longDat, aes(x = cutDat, y = residuals)) +
  geom_boxplot() +
  theme_cowplot(18) +
  geom_hline(yintercept = 0, color = "red" ,linetype = "dashed") +
  xlab("Read depth ") +
  ylab("RRBS DNAm - array DNAm") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


save(ResidualsRRBSArrayPlot, arrayRRBSDistribPlot, arrayRRBSDistribDat,
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/arrayRRBSDistribution.Rdata")


### Supplementary Figure 6 and 8 #####################
### SupplementaryAccuracyVSTime 
load("/mnt/data1/Thea/Simulations/data/fromISCA/benchmarkingDataAll.Rdata")
library(ggplot2)
library(cowplot)
library(reshape2)

colnames(benchmarkResults) = c(1000, seq(2500, 50000, 2500))

nPropDat = melt(benchmarkResults[1:50,])
rDat = melt(benchmarkResults[51:100,])
colnames(rDat) = colnames(nPropDat) = c("var1", "nSites", "value")

nPropPlot = ggplot(nPropDat, aes(x = nSites, y = value, group = nSites)) +
  geom_boxplot() +
  theme_cowplot(18) +
  xlab("Number of sites used to calculate proportion") +
  ylab("Proportion") +
  geom_hline(yintercept = nProp, col = "red", linetype = "dashed")

rPropPlot = ggplot(rDat, aes(x = nSites, y = value, group = nSites)) +
  geom_boxplot() +
  theme_cowplot(18) +
  xlab("Number of sites used to calculate r") +
  ylab("r parameter") +
  geom_hline(yintercept = r, col = "red", linetype = "dashed")

save(r, nProp, nPropPlot,rPropPlot, 
     file = "/mnt/data1/Thea/Simulations/data/dataForPlots/SupplementaryAccuracyVSTime.Rdata")

### Supplementary Figure 7 ############################
### suppOptimisePrior 
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
DNAm = allOrdered[,2:126]

## calculate proportions
priorCalc = function(RRBSm){
  sumCalc = function(RRBSmRow){
    g1 = sum(RRBSmRow < 5, na.rm = T)
    g3 = sum(RRBSmRow > 95, na.rm = T)
    tot = sum(!is.na(RRBSmRow))
    return(c(g1, g3, tot))
  }
  
  propDat = colSums(t(apply(RRBSm, 1, sumCalc)))
  
  priors = c(propDat[1]/propDat[3],
             (propDat[3]-(propDat[1]+propDat[2]))/propDat[3],
             propDat[2]/propDat[3])
  
  return(priors)
}

repfunc = function(nSub){
  x = replicate(n = 50, priorCalc(DNAm[sample(nrow(DNAm), nSub),]))
  return(x)
}


a = sapply(c(1000, seq(5000, 50000, 5000)), repfunc)

trueVal = priorCalc(DNAm)

save(a, trueVal, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/suppOptimisePrior.Rdata")

load("/mnt/data1/Thea/Simulations/data/dataForPlots/suppOptimisePrior.Rdata")
colnames(a) = c(1000, seq(5000, 50000, 5000))
library(reshape2)
library(colortools)

cols = triadic("#A40061")

prior1 = melt(a[seq(1,150,3),])
prior2 = melt(a[seq(2,150,3),])
prior3 = melt(a[seq(3,150,3),])

prior1[,1] = "1"
prior2[,1] = "2"
prior3[,1] = "3"

colnames(prior1) = colnames(prior2) = colnames(prior3) =  c("priorGroup", "subset", "value")

p1 = ggplot(prior1, aes(group = subset, y = value, x = subset)) +
  geom_boxplot(fill = cols[1]) +
  theme_cowplot(18) +
  labs(x = "Number of sites", y = "Prior") +
  geom_hline(yintercept = trueVal[1], linetype = "dashed", col = "red")
p2 = ggplot(prior2, aes(group = subset, y = value, x = subset)) +
  geom_boxplot(fill = cols[2]) +
  theme_cowplot(18) +
  labs(x = "Number of sites", y = "Prior") +
  geom_hline(yintercept = trueVal[2], linetype = "dashed", col = "red")
p3 = ggplot(prior3, aes(group = subset, y = value, x = subset)) +
  geom_boxplot(fill = cols[3]) +
  theme_cowplot(18) +
  labs(x = "Number of sites", y = "Prior") +
  geom_hline(yintercept = trueVal[3], linetype = "dashed", col = "red")

## make legend for across all plots
pLeg = rbind.data.frame(prior1, prior2, prior3)
plotForLeg =
  ggplot(pLeg, aes(x = as.factor(subset), y = value, fill = as.factor(priorGroup))) +
  geom_boxplot() +
  scale_fill_manual(values = cols, labels = c("x<0.05", "0.05<x<0.95", "x>0.95")) +
  theme_cowplot(18) +
  labs(x = "Number of sites", y = "Prior", fill = "DNAm bin")

leg = get_legend(plotForLeg)

save(leg, p1,p2,p3, file = "/mnt/data1/Thea/Simulations/data/dataForPlots/suppOptimisePriorPlot.Rdata")



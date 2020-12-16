#######################################################
### Scripts used to saving plot data in RRBS paper ####
### Dorothea Seiler Vellame ###########################
### ds420@exeter.ac.uk ################################
#######################################################
library(cowplot)
library(ggplot2)

### Figure 1 ##########################################
## RRBSDistribution
load("/mnt/data1/Thea/Simulations/data/dataForPlots/RRBSDistribution.Rdata")
x = plotCovRRBS200$data

cov = ggplot(data = x, aes(x = value)) +
  geom_line(aes(group = X2, y = ..count../1000000), stat="density", alpha=0.2) +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  xlab("Read depth") +
  ylab("Number of DNAm points (million)")

x = plotMethRRBS$data
meth = ggplot(data = x, aes(x = value)) +
  geom_line(aes(group = Var2, y = ..count../1000000), stat="density", alpha=0.2) +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  xlab("DNAm") +
  ylab("Number of DNAm points (million)")


pdf("/mnt/data1/Thea/Simulations/simulationImages/combo/Figure1.pdf", height = 6, width = 12)
plot_grid(cov, meth, labels = "AUTO", nrow = 1, ncol = 2, axis = "b", align = "h")
dev.off()



### Figure 2 ##########################################
load("/mnt/data1/Thea/Simulations/data/dataForPlots/arrayVSRRBSWithFiltering.Rdata")

## make filtered plot
x = RRBSArrayCorPlot$data

## find peak of pearson rd filter
y = RRBSArrayCorLinePlot$data[which(RRBSArrayCorLinePlot$data$type == "Pearson"),]
peakFilt = y[which(y[,2] == max(y[,2])),1]

plotDat = x[which(x$Coverage > peakFilt),]

nonFilt = ggplot(x, aes(x = RRBSmethlong, y = arraylong, col = Coverage)) +
  geom_point() +
  labs(colour = "Read depth") +
  theme_cowplot(18) +
  xlab("RRBS DNAm") +
  ylab("Array DNAm") +
  scale_color_gradient(low="blue", high="red", limits = c(0,125), name = "Read depth")
filt = ggplot(plotDat, aes(x = RRBSmethlong, y = arraylong, col = Coverage)) +
  geom_point() +
  labs(colour = "Read depth") +
  theme_cowplot(18) +
  xlab("RRBS DNAm") +
  ylab("Array DNAm") +
  scale_color_gradient(low="blue", high="red", limits = c(0,125), name = "Read depth")

noLeg = plot_grid(nonFilt + theme(legend.position="none"),
                  filt + theme(legend.position="none"),
                  labels = c("D","E"), ncol = 2, nrow = 1)
legend <- get_legend(
  nonFilt + theme(legend.box.margin = margin(0, 0, 0, 12))
)

corPlots = plot_grid(noLeg, legend, rel_widths = c(2, .3))


load("/mnt/data1/Thea/Simulations/data/rdCorrDNAmRealSimPlot.Rdata")

p1 = ggplot(data = RRBSArrayCorLinePlot$data, aes(x = covVals, y = corr, col = type)) +
  geom_line(size = 2, col = "#F8766D") +
  theme_cowplot(18) +
  xlab("Read depth filter") +
  ylab("Correlation")  +
  theme(legend.title = element_blank())
p2 = ggplot(data = RRBSArrayRMSEPlot$dat, aes(x = covVals, y = corr, col = type)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#1A48C4")) +
  theme_cowplot(18) +
  xlab("Read depth filter") +
  ylab("RMSE") +
  theme(legend.title = element_blank())

linePlots = 
  plot_grid(rdSimCorPlot,
            rdSimErrorPlot +ylab("RMSE"),
            p1,
            p2 + theme(legend.position="none"),
            labels = c("Ci", "Cii", "Fi", "Fii"), nrow = 4, ncol = 1, axis = "lr", align = "v")

## load violin plots and rd with extremes
load("/mnt/data1/Thea/Simulations/data/dataForPlots/extremesPlot.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/simulatedPrecision.Rdata")
topright = plot_grid(extremesPlot, simulatedPrecision, ncol = 2, labels = "AUTO", rel_widths = c(1,1.3))
rightPlots = plot_grid(topright, corPlots, nrow = 2)

pdf("/mnt/data1/Thea/Simulations/simulationImages/combo/Figure2.pdf", height = 11, width = 20)
plot_grid(rightPlots, linePlots, nrow = 1, rel_widths = c(1, 0.35))
dev.off()





### Figure 3 ##########################################
load("/mnt/data1/Thea/Simulations/data/dataForPlots/OverlapBetweenSampleCpGs.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/corrHistBetweenSamples.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/dataRemaining.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/rdCorCol.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/annotationPlot.Rdata")

wideplot = plot_grid(propRemainPlot +ylab("Proportion remaining"),
                     overlapBetweenSamplesPlot + ylim(c(0,0.65)) +ylab("Proportion present")+
                       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
                     ncol = 1,
                     axis = "bl",
                     align = "hv",
                     labels = c("A", "C"))
otherPlot = plot_grid(rdCorHist,
                      rdCorColPlot + labs(col = "DNAm\npoints\nper site"),
                      labels = c("B", "D"),
                      ncol = 1,
                      axis = "blr",
                      align = "hv")

pdf("/mnt/data1/Thea/Simulations/simulationImages/combo/Figure3.pdf", height = 9, width = 20)
plot_grid(wideplot, otherPlot, annotationPlot,
          labels = c("","","E"),
          rel_widths = c(1.5,1,0.8),
          ncol = 3,
          axis = "bl",
          align = "hv")
dev.off()


### Figure 4 ##########################################
load("/mnt/data1/Thea/Simulations/data/dataForPlots/simPlotData.Rdata")

rdPlot =
  ggplot(rdPlot$data, aes(x = rd, y = dat, linetype = group)) +
  geom_line(size = 1.5, col = "#1AC496") +
  ylim(0,100) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  xlab("Read depth") +
  ylab("Power (%)") +
  theme_cowplot(18) +
  scale_linetype_discrete(name = "ΔμDNAm", labels = c("0.2", "0.1", "0.05"))


mdPlot = ggplot(mdPlot$data, aes(x = md, y = dat, linetype = group)) +
  geom_line(size = 1.5, col = "#1A48C4") +
  ylim(0,100) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  xlab("Mean DNAm difference") +
  ylab("Power (%)") +
  theme_cowplot(18) +
  scale_linetype_discrete(name = "N per\ngroup", labels = c("20", "10", "5"))

muPlot = ggplot(muPlot$data, aes(x = mu/100, y = dat, linetype = group)) +
  geom_line(size = 1.5, col = "#9DC41A") +
  ylim(0,100) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  theme_grey(base_size = 18) +
  xlab("Mean DNAm") +
  ylab("Power (%)") +
  theme_cowplot(18) +
  scale_linetype_discrete(name = "μRD", labels = c("50", "30", "10"))

varPlot = ggplot(muPlot$data, aes(x = mu/100, y = var, linetype = group)) +
  geom_line(size = 1.5, col = "#9DC41A") +
  theme_grey(base_size = 18) +
  xlab("Mean DNAm") +
  ylab("Variance") +
  theme_cowplot(18) +
  scale_linetype_discrete(name = "μRD", labels = c("50", "30", "10")) 

sampleSizePlot =
  ggplot(sampleSizePlot$data, aes(x = meanNPoss, y = power, col = balanceGroup)) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  ylim(0,100) +
  xlab("Total sample size") +
  ylab("Power (%)") +
  labs(col = expression(atop("Sample", paste("split")))) +
  theme_cowplot(18)

topPlots = plot_grid(NULL, rdPlot, sampleSizePlot, NULL, labels = c("", "A", "B", ""),
                     rel_widths = c(0.5,1,1.1,0.5), ncol = 4, align = "h")
bottomPlots = plot_grid(mdPlot, muPlot, varPlot, labels = c("C","D","E"),
                        ncol = 3, axis = "b", align = "h")

grDevices::cairo_pdf("/mnt/data1/Thea/Simulations/simulationImages/combo/Figure4.pdf", height = 8, width = 15)
plot_grid(topPlots, bottomPlots, nrow = 2)
dev.off()


### Figure 5 ##########################################
## load packages
library(ggplot2)
library(cowplot)
library(reshape2)

load("/mnt/data1/Thea/Simulations/data/dataForPlots/simQQ.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/powerWithPOWEREDBiSeq.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/diffInPower1_75rd.Rdata")
load("/mnt/data1/Thea/Simulations/data/dataForPlots/minRDPlot.Rdata")

qqplots = plot_grid(RDQQFull, minRDPlot, DNAmQQ, ncol = 3,
                    labels = c("Ai", "Aii", "B"),
                    axis = "l", align = "tb")

rightplots = plot_grid(powerWithPOWEREDBiSeq + labs(col = "Samples\nneeded", shape = "Samples\nneeded"),
                       powerDiffPlot + labs(x = "Samples needed", y = "Difference in\npower (%)"), ncol = 2,
                       labels = c("C", "D"),
                       axis = "tb", align = "h", rel_widths = c(1, 0.5))



pdf("/mnt/data1/Thea/Simulations/simulationImages/combo/Figure5.pdf", height = 7, width = 11)
plot_grid(qqplots, rightplots, ncol = 1, rel_heights = c(1, 1))
dev.off()



### Supplementary Figure 1 ############################

load("/mnt/data1/Thea/Simulations/data/dataForPlots/arrayRRBSDistribution.Rdata")
pdf("/mnt/data1/Thea/Simulations/simulationImages/supplementary/arrayRRBSDistrib.pdf",width = 9, height = 7)
ggplot(data = arrayRRBSDistribDat, aes(x = meth, col = type)) +
  geom_line(aes(group = sample, y = ..count../1000), stat="density", alpha = 0.3) +
  theme_cowplot(18) +
  theme(legend.title = element_blank()) +
  xlab("DNAm") +
  ylab("Number of sites (thousand)") +
  scale_color_manual(values=c('#2a18a2','#90a218')) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) 
dev.off()

### Supplementary Figure 2 ############################
load("/mnt/data1/Thea/Simulations/data/dataForPlots/arrayRRBSDistribution.Rdata")
pdf("/mnt/data1/Thea/Simulations/simulationImages/supplementary/residualsRRBSArray.pdf",width = 10, height = 7)
ResidualsRRBSArrayPlot
dev.off()

### Supplementary Figure 5 ############################
## power for meanDiff 0.06 nSN = 60, rd = 15 
library(ggplot2)
library(cowplot)
library(reshape2)
library(tidyr)

## see Figure 5 for the making of powerOutHistDat

y = matrix(ncol = 9, nrow = 0)
for (i in 0:9){
  load(paste("/mnt/data1/Thea/Simulations/data/fromISCA/powerOutHistDat_", i,".Rdata", sep = ""))
  x = melt(x)
  y = rbind(y, spread(x, key = "Var1", value = "value"))
}

pdf("/mnt/data1/Thea/Simulations/simulationImages/supplementary/PBSoutputs.pdf", height = 6, width = 6)
ggplot(y, aes(x = power)) +
  geom_histogram(bins = 50) +
  theme_cowplot(18) +
  labs(x = "Power (%)", y = "Number of POWEREDBiSeq\ncalculations")
dev.off()


### Supplementary Figure 6 and 8 ############################
load("/mnt/data1/Thea/Simulations/data/dataForPlots/SupplementaryAccuracyVSTime.Rdata")
library(cowplot)
library(ggplot2)
pdf("/mnt/data1/Thea/Simulations/simulationImages/supplementary/rAccuractBoxplots.pdf", height = 7, width = 9)
rPropPlot + labs(y = "r", x = "Number of sites")
dev.off()
pdf("/mnt/data1/Thea/Simulations/simulationImages/supplementary/nAccuractBoxplots.pdf", height = 7, width = 9)
nPropPlot + labs(x = "Number of sites")
dev.off()

### Supplementary Figure 7 ############################
load("/mnt/data1/Thea/Simulations/data/dataForPlots/suppOptimisePriorPlot.Rdata")
priorPlots = plot_grid(p1,p2,p3, ncol = 1)

pdf("/mnt/data1/Thea/Simulations/simulationImages/supplementary/suppOptimisePrior.pdf", height = 10, width = 12)
plot_grid(priorPlots, leg, ncol = 2, rel_widths = c(1,0.2))
dev.off()




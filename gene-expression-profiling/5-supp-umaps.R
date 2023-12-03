###  
### Supplementary UMAPs
###

# Set random seed
set.seed(83539)

# Set margins
mmm = c(12,2,2,10)

# Load normalised TCGA data and annotation
tcga_normalised = as.matrix(read.csv("gene-expression-profiling/tcga-normalised.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))
molecular_data = readxl::read_excel("data/Supplementary Tables.xlsx", sheet=2)
molecular_data = molecular_data[match(colnames(tcga_normalised),molecular_data$`Patient ID`),]

# Load library
library(umap)

# Select group
w = which(molecular_data$GEP == "Class I")

# Perform UMAP
u = umap(t(tcga_normalised[,w]), preserve.seed = T)

# Filename
filename = "gene-expression-profiling/supp-umap-1.png"

# Generate file
file = paste0(filename)
png(file, res=600, 3000, 2700)  
par(mar=mmm)
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)

# Determine colors
colors = c("#DDDDDD","#777777")
values = rep(colors[1], length(w))
new = molecular_data[w,]
values[which(!is.na(new$`EIF1AX mutation`))] = colors[2]

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

yat = c(-3,-2,-1,0,1,2)
xat = c(-4,-3,-2,-1,0,1,2,3,4,5)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout, pch=16,col=values)

line = (ylim[2]-ylim[1]) / 9
text(mean(xlim), ylim[2], cex=1.1,"GEP class I only (n=40)", xpd=T,pos=3)

x = xlim[2] + (xlim[2] - xlim[1]) /64
y = mean(ylim) + 1 * line
text(x, y, cex=1.1,pos=4,"EIF1AX", font=3, xpd=T)

x = xlim[2] + (xlim[2] - xlim[1]) / 16
y = mean(ylim)
points(x,y, pch=16,col=colors[1],xpd=T)
text(x, y, cex=1.1,pos=4,"wild-type", xpd=T)

y = mean(ylim) - 1 * line
points(x,y, pch=16,col=colors[2],xpd=T)
text(x, y, cex=1.1,pos=4,"mutant", xpd=T)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

dev.off()
system(paste("open", file))

###

# Filename
filename = "gene-expression-profiling/supp-umap-2.png"

# Generate file
file = paste0(filename)
png(file, res=600, 3000, 2700)  
par(mar=mmm)
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)

# Determine colors
colors = c("#DDDDDD","#777777")
values = rep(colors[1], length(w))
new = molecular_data[w,]
values[which(!is.na(new$`SF3B1 mutation`))] = colors[2]

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

yat = c(-3,-2,-1,0,1,2)
xat = c(-4,-3,-2,-1,0,1,2,3,4,5)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout, pch=16,col=values)

line = (ylim[2]-ylim[1]) / 9
text(mean(xlim), ylim[2], cex=1.1,"GEP class I only (n=40)", xpd=T,pos=3)

x = xlim[2] + (xlim[2] - xlim[1]) /64
y = mean(ylim) + 1 * line
text(x, y, cex=1.1,pos=4,"SF3B1", font=3, xpd=T)

x = xlim[2] + (xlim[2] - xlim[1]) / 16
y = mean(ylim)
points(x,y, pch=16,col=colors[1],xpd=T)
text(x, y, cex=1.1,pos=4,"wild-type", xpd=T)

y = mean(ylim) - 1 * line
points(x,y, pch=16,col=colors[2],xpd=T)
text(x, y, cex=1.1,pos=4,"mutant", xpd=T)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

dev.off()
system(paste("open", file))

###

# Filename
filename = "gene-expression-profiling/supp-umap-3.png"

# Generate file
file = paste0(filename)
png(file, res=600, 3000, 2700)  
par(mar=mmm)
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)

# Determine colors
colors = c("#DDDDDD","#777777","#333333")
values = rep(colors[1], length(w))
new = molecular_data[w,]
values[which(new$`Chr. 8q (DNA)` == "gain")] = colors[2]
values[which(stringr::str_detect(new$`Chr. 8q (DNA)`,"amplification"))] = colors[2]

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

yat = c(-3,-2,-1,0,1,2)
xat = c(-4,-3,-2,-1,0,1,2,3,4,5)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout, pch=16,col=values)

line = (ylim[2]-ylim[1]) / 9
text(mean(xlim), ylim[2], cex=1.1,"GEP class I only (n=40)", xpd=T,pos=3)

x = xlim[2] + (xlim[2] - xlim[1]) /64
y = mean(ylim) + 1 * line
text(x, y, cex=1.1,pos=4,"Chromosome 8q", font=1, xpd=T)

x = xlim[2] + (xlim[2] - xlim[1]) / 16
y = mean(ylim)
points(x,y, pch=16,col=colors[1],xpd=T)
text(x, y, cex=1.1,pos=4,"disomy", xpd=T)

y = mean(ylim) - 1 * line
points(x,y, pch=16,col=colors[2],xpd=T)
text(x, y, cex=1.1,pos=4,"gain/amplification", xpd=T)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

dev.off()
system(paste("open", file))

###
###
###


# Select group
w = which(molecular_data$GEP == "Class II")

# Perform UMAP
u = umap(t(tcga_normalised[,w]), preserve.seed = T)

# Filename
filename = "gene-expression-profiling/supp-umap-4.png"

# Generate file
file = paste0(filename)
png(file, res=600, 3000, 2700)  
par(mar=mmm)
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)

# Determine colors
colors = c("#DDDDDD","#777777")
values = rep(colors[1], length(w))
new = molecular_data[w,]
values[which(!is.na(new$`EIF1AX mutation`))] = colors[2]

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

yat = c(-3,-2,-1,0,1,2)
xat = c(-4,-3,-2,-1,0,1,2,3,4,5)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout, pch=16,col=values)

line = (ylim[2]-ylim[1]) / 9
text(mean(xlim), ylim[2], cex=1.1,"GEP class II only (n=40)", xpd=T,pos=3)

x = xlim[2] + (xlim[2] - xlim[1]) /64
y = mean(ylim) + 1 * line
text(x, y, cex=1.1,pos=4,"EIF1AX", font=3, xpd=T)

x = xlim[2] + (xlim[2] - xlim[1]) / 16
y = mean(ylim)
points(x,y, pch=16,col=colors[1],xpd=T)
text(x, y, cex=1.1,pos=4,"wild-type", xpd=T)

y = mean(ylim) - 1 * line
points(x,y, pch=16,col=colors[2],xpd=T)
text(x, y, cex=1.1,pos=4,"mutant", xpd=T)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

dev.off()
system(paste("open", file))

###

# Filename
filename = "gene-expression-profiling/supp-umap-5.png"

# Generate file
file = paste0(filename)
png(file, res=600, 3000, 2700)  
par(mar=mmm)
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)

# Determine colors
colors = c("#DDDDDD","#777777")
values = rep(colors[1], length(w))
new = molecular_data[w,]
values[which(!is.na(new$`SF3B1 mutation`))] = colors[2]

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

yat = c(-3,-2,-1,0,1,2)
xat = c(-4,-3,-2,-1,0,1,2,3,4,5)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout, pch=16,col=values)

line = (ylim[2]-ylim[1]) / 9
text(mean(xlim), ylim[2], cex=1.1,"GEP class II only (n=40)", xpd=T,pos=3)

x = xlim[2] + (xlim[2] - xlim[1]) /64
y = mean(ylim) + 1 * line
text(x, y, cex=1.1,pos=4,"SF3B1", font=3, xpd=T)

x = xlim[2] + (xlim[2] - xlim[1]) / 16
y = mean(ylim)
points(x,y, pch=16,col=colors[1],xpd=T)
text(x, y, cex=1.1,pos=4,"wild-type", xpd=T)

y = mean(ylim) - 1 * line
points(x,y, pch=16,col=colors[2],xpd=T)
text(x, y, cex=1.1,pos=4,"mutant", xpd=T)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

dev.off()
system(paste("open", file))

###

# Filename
filename = "gene-expression-profiling/supp-umap-6.png"

# Generate file
file = paste0(filename)
png(file, res=600, 3000, 2700)  
par(mar=mmm)
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)

# Determine colors
colors = c("#DDDDDD","#777777","#333333")
values = rep(colors[1], length(w))
new = molecular_data[w,]
values[which(new$`Chr. 8q (DNA)` == "gain")] = colors[2]
values[which(stringr::str_detect(new$`Chr. 8q (DNA)`,"amplification"))] = colors[2]

plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")

yat = c(-3,-2,-1,0,1,2)
xat = c(-4,-3,-2,-1,0,1,2,3,4,5)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout, pch=16,col=values)

line = (ylim[2]-ylim[1]) / 9
text(mean(xlim), ylim[2], cex=1.1,"GEP class II only (n=40)", xpd=T,pos=3)

x = xlim[2] + (xlim[2] - xlim[1]) /64
y = mean(ylim) + 1 * line
text(x, y, cex=1.1,pos=4,"Chromosome 8q", font=1, xpd=T)

x = xlim[2] + (xlim[2] - xlim[1]) / 16
y = mean(ylim)
points(x,y, pch=16,col=colors[1],xpd=T)
text(x, y, cex=1.1,pos=4,"disomy", xpd=T)

y = mean(ylim) - 1 * line
points(x,y, pch=16,col=colors[2],xpd=T)
text(x, y, cex=1.1,pos=4,"gain/amplification", xpd=T)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

dev.off()
system(paste("open", file))


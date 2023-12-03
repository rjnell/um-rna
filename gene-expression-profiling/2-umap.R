###  
### STEP 2 - VISUALISE UMAP
###

# Set random seed
set.seed(83539)

# Load normalised TCGA data and annotation
tcga_normalised = as.matrix(read.csv("gene-expression-profiling/tcga-normalised.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))
molecular_data = readxl::read_excel("data/Supplementary Tables.xlsx", sheet=2)
molecular_data = molecular_data[match(colnames(tcga_normalised),molecular_data$`Patient ID`),]

# Load library
library(umap)

# Generate UMAP
u = umap(t(tcga_normalised), preserve.seed = T)
plot(u$layout, col=RColorBrewer::brewer.pal(9,"Paired")[as.factor(molecular_data$GEP)], pch=16)

# Generate figure
file = paste0("gene-expression-profiling/figure-1A.png")
png(file, res=600, 3000, 2700)  
par(mar=c(7,6,5,4))
xlim = c(min(u$layout[,1])-(max(u$layout[,1])-min(u$layout[,1]))*0.25,max(u$layout[,1])+(max(u$layout[,1])-min(u$layout[,1]))*0.25)
ylim = c(min(u$layout[,2])-(max(u$layout[,2])-min(u$layout[,2]))*0.25,max(u$layout[,2])+(max(u$layout[,2])-min(u$layout[,2]))*0.25)
ccolors = c("#AAAAAA","#777777")
names(ccolors) = c("Class I","Class II")
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
points(u$layout, pch=16,col=ccolors[molecular_data$GEP])
yy=ylim[1]-(ylim[2]-ylim[1])/16*5
points(xlim[1]+(xlim[2]-xlim[1])/16*1,yy, pch=16,col=ccolors["Class I"],xpd=T)
text(xlim[1]+(xlim[2]-xlim[1])/16*1,yy, cex=1.1,pos=4,"GEP class I", xpd=T)
points(xlim[1]+(xlim[2]-xlim[1])/16*9,yy, pch=16,col=ccolors["Class II"],xpd=T)
text(xlim[1]+(xlim[2]-xlim[1])/16*9,yy, cex=1.1,pos=4,"GEP class II", xpd=T)
segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)
dev.off()
system(paste("open", file))
###
##  Figure 4A
###

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
lumc_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=3)


colors = list()
colors[["female"]] = RColorBrewer::brewer.pal(12, "Paired")[5]
colors[["Class I"]] = "#AAAAAA"
colors[["Class II"]] = "#777777"
colors[["PLCB4"]] = RColorBrewer::brewer.pal(8,"YlGnBu")[4]
colors[["CYSLTR2"]] = RColorBrewer::brewer.pal(8,"YlGnBu")[6]
colors[["EIF1AX-exon1"]] = RColorBrewer::brewer.pal(10, "Paired")[10]
colors[["EIF1AX-exon2"]] = RColorBrewer::brewer.pal(9, "Paired")[9]
colors[["SF3B1-hotspot"]] = RColorBrewer::brewer.pal(10, "Paired")[8]
colors[["SF3B1-other"]] = RColorBrewer::brewer.pal(10, "Paired")[7]
colors[["BAP1"]] = RColorBrewer::brewer.pal(12, "Paired")[12]

colors[["loss"]] = RColorBrewer::brewer.pal(12, "Paired")[6]
colors[["gain"]] = RColorBrewer::brewer.pal(11, "RdYlGn")[9]
colors[["amp"]] = RColorBrewer::brewer.pal(11, "RdYlGn")[10]

colors[["DNA + RNA"]] = "#999999"
colors[["RNA"]] = "#2196F3"
colors[["DNA"]] = "#F44336"

m = 82
bar = "darkblue"
s = 87
z = 10

i_stop = 90

file = "lumc-cohort/figure-4A.png"
png(file, res=600, width=5300, height=7500)
par(mar=c(5,5,5,5))
plot(c(-15,95), 
     c(-10,40),
     type="n", 
     axes=F, 
     xaxs="i", 
     yaxs="i", 
     xlab="", 
     ylab="") 

# Set seed
set.seed(1234)

# Load lumc data
#lumc_normalised = as.matrix(read.csv("../um-heterogeneity/res/lumc_normalized_80.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))

# Load classifier data
classifier = read.csv("gene-expression-profiling/classifier.tsv", stringsAsFactors = F, sep = "\t")
# Set random seed
set.seed(83539)
# Plot 
lumc_classifier = lumc_normalised[classifier$ENSG,lumc_data$Tumor_ID]
lumc_classifier = lumc_classifier[which(matrixStats::rowSds(lumc_classifier)>0.5,),]
lumc_classifier_ph = pheatmap::pheatmap(lumc_classifier, 
                                        scale="row", 
                                        custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                                        cutree_cols = 2, 
                                        cutree_rows = 2,
                                        silent = T)
ord = rev(lumc_classifier_ph$tree_col$labels[lumc_classifier_ph$tree_col$order])
lumc_data = lumc_data[match(ord, lumc_data$Tumor_ID),]
lumc_classifier = lumc_classifier[lumc_classifier_ph$tree_row$order,match(ord, colnames(lumc_classifier))]


rows = 1:nrow(lumc_classifier)
pos=c(s,s+2.5,s+5,s+7.5, s+10)
#segments(pos, 29+0.4, pos, 29-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#eeeeee", xpd=T)
segments(pos, 14+1.4, pos, 14-0.4, lwd=1.4, col="#eeeeee", xpd=T)

i_stop = 3
y = 33.4
xi= 0.4
yi = 6.8/length(rows)
palette = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(101))
for (row in rows) {
  for (i in 1:5) {
    x = i+1*(i-1) 
    if (i > i_stop) { x = x + 2 }
    col = round((lumc_classifier[row, i] - min(lumc_classifier[row, ])) / (max(lumc_classifier[row, ]) - min(lumc_classifier[row, ])) * 100) + 1
    rect(x-xi, y-yi/2, x+xi+1, y+yi/2, border=NA, col=palette[col],xpd=T)
  }
  y = y - yi
}







white = "#EEEEEE"
border = NA
xi = 0.4

y = 33-7
for (i in 1:5) {
  col = white
  x = i+1*(i-1) 
  if (i > i_stop) { x = x + 2 }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=colors[[lumc_data$Inferred_GEP[i]]])
}


y = 30

y=25
for (gaq in c("GNAQ", "GNA11", "CYSLTR2", "PLCB4")) {
  # SF3B1
  y = y-1
  for (i in 1:5) {
    x = i+1*(i-1) 
    if (i > i_stop) { x = x + 2 }
    val = lumc_data$Gaq_mutation[i]
    col = white
    if (!is.na(val) & stringr::str_detect(val, gaq)) {
      col = "#777777"
    }
    rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  }
}

# EIF1AX
y = 17
for (i in 1:5) {
  x = i+1*(i-1) 
  if (i > i_stop) { x = x + 2 }
  val = lumc_data$EIF1AX_mutation[i]
  col = white
  if (!is.na(val)) {
    col = "#777777"
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
}

# SF3B1
y = 18
for (i in 1:5) {
  x = i+1*(i-1) 
  if (i > i_stop) { x = x + 2 }
  val = lumc_data$SF3B1_mutation[i]
  col = white
  if (!is.na(val)) {
    col = "#777777"
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
}

# BAP1
y = 19
for (i in 1:5) {
  x = i+1*(i-1) 
  if (i > i_stop) { x = x + 2 }
  val = lumc_data$BAP1_mutation[i]
  col = white
  if (!is.na(val)) {
    col = "#777777"
  }
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
}

# Chr3p
y = 15
for (i in 1:5) {
  x = i+1*(i-1) 
  if (i > i_stop) { x = x + 2 }
  val = lumc_data$Chr3p_copy_number[i]
  col = white
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  if (!is.na(val)) {
    col = white
    if (val == "trisomy 3") { 
      col = "#AAAAAA"
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    } 
    else {
      if (val == "loss") { 
        col = "#777777"
        rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
      }
      else if (val == "isodisomy 3") { 
        col = "#777777"
        rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
      }
      else if (val == "partial monosomy 3") { 
        col = "#F44336"
        #Arect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white)
        rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col="#777777")
      } 
      
    }
  }
}

# Chr8q
y = 14
for (i in 1:5) {
  x = i+1*(i-1) 
  if (i > i_stop) { x = x + 2 }
  val = lumc_data$Chr8q_copy_number[i]
  col = white
  rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
  if (!is.na(val)) {
    col = white
    if (val == "gain") { 
      col = "#AAAAAA"
    }
    else if (val == "amplification") { 
      col = "#777777"
    }
    else if (val == "amplification (biallelic)") { 
      col = "#777777"
    } 
    rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=col)
    
    text(x,y-1,cex=0.9,srt=90,adj=1,labels=lumc_data$ID[i], col="#333333")
  }
}





pos = 1-0.7-2

limits=c(33,26)
axis(side = 2, at = mean(limits), labels = "GEP classifier", las=2, cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
segments(pos, limits[1]+0.4, pos, limits[2]-0.4, lwd=1.4, col="#b1b1b1", xpd=T)

axis(side = 2, at = 24:21, labels = c("GNAQ", "GNA11", "CYSLTR2", "PLCB4"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)



axis(side = 2, at = 17:19, labels = c("EIF1AX","SF3B1","BAP1"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
axis(side = 2, at = 15:14, labels = c("Chromosome 3","Chromosome 8q"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=1, col.ticks = "#b1b1b1", lwd.ticks = 1.4)

#pos = c(pos, s)
segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 14+1.4, pos, 14-0.4, lwd=1.4, col="#b1b1b1")


xn = -22.5+22.5
xn = 20
y = 7+4
y = 33
x =xn+1
text(x-2.5, y, labels = "GEP classifier", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#AAAAAA", xpd=T)
text(x+1, y, labels = "class I", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#777777", xpd=T)
text(x+1, y, labels = "class II", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

y=y+1
x=x+17.55
text(x-2.3, y, labels = "low", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
text(x+18, y, labels = "high", cex=1.1, col="#333333", font=1, pos=2, xpd=T)

palette = colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(101)
for (k in 1:101) {
  rect(x+5+4.5/101*(k-1), y-0.2, x+5+4.5/101*(k), y+0.2, border=NA, col=palette[k],xpd=T)
}
rect(x+5, y-0.2, x+9.5, y+0.2, border=border, xpd=T)
text(x+7.5, y-1, labels = "expression", cex=1.1, col="#333333", font=1, xpd=T)

y=y-7
x=xn+1
text(x-2.5, y, labels = "Mutations", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=white, xpd=T)
text(x+1, y, labels = "wild-type (RNA + DNA)", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn+1
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#777777", xpd=T)
text(x+1, y, labels = "mutant (RNA + DNA)", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


#y = yn
#y = 7+4
#xn= 82-39
y = y-2
x=xn
text(x-2.5, y, labels = "Chromosome 3 (RNA / DNA)", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=white, xpd=T)
text(x+1, y, labels = "no LOH / no alteration", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#x=xn+15
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#777777", xpd=T)
text(x+1, y, labels = "LOH / loss", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#x=xn+15


y=y-2
text(x-2.5, y, labels = "Chromosome 8q (RNA / DNA)", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col=white, xpd=T)
text(x+1, y, labels = "no LOH / no alteration", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#777777", xpd=T)
text(x+1, y, labels = "LOH / gain", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#x=xn+15
y = y-1
rect(x-xi, y-0.4, x+xi+1, y+0.4, border=border, col="#aaaaaa", xpd=T)
text(x+1, y, labels = "LOH / amplification", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
#x=xn+15

dev.off()
system(paste("open", file))





# Load dependencies (pheatmap)
source("bin/pheatmap.r")
library(grid)
library(RColorBrewer)
library(scales)
library(gtable)
library(stats)
library(grDevices)
library(graphics)

# Set random seed
set.seed(83539)

# Plot and save dendrogram
file = "lumc-cohort/figure-4A-dendrogram.png"
png(file, res=600, width=1300, height=1500)
dendrogram = pheatmap(lumc_classifier, 
                      scale="row", 
                      custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                      cutree_cols = 2, 
                      show_rownames = F,
                      treeheight_col=20)
dev.off()
system(paste("open", file))


lumc_classifier_ph

library(umap)
# Set random seed
set.seed(83539)

# Generate UMAP
u = umap(t(lumc_normalised[,lumc_data$Tumor_ID]), preserve.seed = T, n_neighbors=4)
plot(u$layout, pch=16)



file = paste0("lumc-cohort/figure-4A-umap.png")
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
#segments(xlim[1],yat,xlim[2],lwd=1.4,col="#EEEEEE", xpd=T)
axis(side = 2, at = ylim[1], labels=c(""),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
axis(side = 1, at = xlim[1], labels=c(""),las=1,col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
mtext(side=1, "UMAP 1",cex = 1.1,col="#333333", line =1)
mtext(side=2, "UMAP 2",cex = 1.1,col="#333333", line =1)
points(u$layout[,1]-0.1, u$layout[,2], pch=16)
text(u$layout[,1]+0.2, u$layout[,2], labels=lumc_data$ID[match(rownames(u$layout),lumc_data$Tumor_ID)])

points(xlim[1]+(xlim[2]-xlim[1])/16*1,-5, pch=16,col=ccolors["Class I"],xpd=T)
text(xlim[1]+(xlim[2]-xlim[1])/16*1,-5, cex=1.1,pos=4,"GEP class I", xpd=T)
points(xlim[1]+(xlim[2]-xlim[1])/16*9,-5, pch=16,col=ccolors["Class II"],xpd=T)
text(xlim[1]+(xlim[2]-xlim[1])/16*9,-5, cex=1.1,pos=4,"GEP class II", xpd=T)



#mtext(side = 2, text = "Reference abundance", line = 5.4, at = mean(yat), cex=1.1)
mtext(side = 2, text = "Digital PCR", line = 3.4, cex=1.1, col="#333333",at=ym/2)
mtext(side = 2, text = "quantification", line = 2.5, cex=1.1, col="#333333",at=ym/2)

segments(xlim[1],ylim[1],xlim[2],ylim[1],lwd=1.4,col="#b1b1b1",xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],lwd=1.4,col="#b1b1b1",xpd=T)

cols = c("#333333","#333333","#333333","#b1b1b1","#b1b1b1", "#333333")#$ RColorBrewer::brewer.pal(9,"Paired")[c(2,4,6)]

k = 1
for (i in c(0.85,0.75,0.35)) {
  cols[k] = rgb(i,i,i)
  k = k +1
}
#col="#808080"

for (i in 1:nrow(data)) {
  j = i*1.5
  rect(j-1.2, as.numeric(data[i,5]), j-0.8, 0, pch=16, cex=0.9, col=cols[1], border=NA, lwd=1.4, xpd=T)
  arrows(j-1, as.numeric(data[i,6]), j-1, as.numeric(data[i,7]), length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  rect(j-0.7, as.numeric(data[i,2]), j-0.3, 0, pch=16, cex=0.9, col=cols[2], border=NA, lwd=1.4, xpd=T)
  arrows(j-0.5, as.numeric(data[i,3]), j-0.5, as.numeric(data[i,4]), length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  
  #text(j-0.5, as.numeric(data[i,4])+ymax/10, labels=paste0(format(as.numeric(data[i,2]), nsmall = 0), ""), cex=0.8, col="#333333", xpd=T)
  
}

#segments(1-0.7, as.numeric(data[1,2]), 2.5, col=col, lwd=1.4, lty=3)
#segments(2-0.7, as.numeric(data[2,2]), 2.5, col=col, lwd=1.4, lty=3)
#arrows(2.5, as.numeric(data[2,2])+ymax/25, 2.5, as.numeric(data[1,2])-ymax/25, length=0.05, code=1, col="#333333", lwd=1.4)


axis(side = 1, at = (0:nrow(data))*1.5, labels=rep("",nrow(data)+1), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, xpd=T, tick = 0)
for (i in 1:nrow(data)) {
  #axis(side = 3, at = i*1.5-0.75, font=2,labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID==data[i,1])], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0, xpd=T)
  text(i*1.5-0.75, ym+1.5,font=2,labels=lumc_data$LUMC_ID[which(lumc_data$Tumor_ID==data[i,1])], col="#333333", cex=1.1, xpd=T)
  #axis(side = 1, at = 1-0.5, labels=expression('var'[1]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
  #axis(side = 1, at = 2-0.5, labels=expression('var'[2]), col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0, tick = 0)
}

#segments(0,0,3*1.5,lwd=1.4,col="#b1b1b1", xpd=T)

dev.off()
system(paste("open", file))




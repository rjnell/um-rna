###
##  FIGURE 1C
###

# Set random seed
set.seed(83539)

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
tcga_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=2)

tcga_data = tcga_data[order(tcga_data$`Gaq detectability`),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "PLCB4")),which(!str_detect(tcga_data$`Gαq mutation`, "PLCB4")),which(is.na(tcga_data$`Gαq mutation`))),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "CYSLTR2")),which(!str_detect(tcga_data$`Gαq mutation`, "CYSLTR2")),which(is.na(tcga_data$`Gαq mutation`))),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "GNA11")),which(!str_detect(tcga_data$`Gαq mutation`, "GNA11")),which(is.na(tcga_data$`Gαq mutation`))),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "GNAQ")),which(!str_detect(tcga_data$`Gαq mutation`, "GNAQ")),which(is.na(tcga_data$`Gαq mutation`))),]

tcga_data = tcga_data[order(tcga_data$`BAP1 detectability`),]
tcga_data = tcga_data[order(tcga_data$`SF3B1 detectability`),]
tcga_data = tcga_data[order(tcga_data$`EIF1AX detectability`),]
tcga_data = tcga_data[order(tcga_data$`Chr. 3 (DNA)`,decreasing = T),]
tcga_data = tcga_data[order(tcga_data$GEP),]


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

file = "gene-expression-profiling/figure-1C.png"
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

# Load TCGA data
tcga_normalised = as.matrix(read.csv("gene-expression-profiling/tcga-normalised.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))

# Load classifier data
classifier = read.csv("gene-expression-profiling/classifier.tsv", stringsAsFactors = F, sep = "\t")

# Plot 
tcga_classifier = tcga_normalised[classifier$ENSG,]
tcga_classifier = tcga_classifier[which(matrixStats::rowSds(tcga_classifier)>0.5,),]
tcga_classifier_ph = pheatmap::pheatmap(tcga_classifier, 
                                        scale="row", 
                                        custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                                        cutree_cols = 2, 
                                        cutree_rows = 2,
                                        silent = T)
ord = rev(tcga_classifier_ph$tree_col$labels[tcga_classifier_ph$tree_col$order])
tcga_data = tcga_data[match(ord, tcga_data$`Patient ID`),]
tcga_classifier = tcga_classifier[tcga_classifier_ph$tree_row$order,match(ord, colnames(tcga_classifier))]


rows = 1:nrow(tcga_classifier)
pos=c(s,s+2.5,s+5,s+7.5, s+10)
#segments(pos, 29+0.4, pos, 29-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#eeeeee", xpd=T)
#segments(pos, 14+1.4, pos, 14-0.4, lwd=1.4, col="#eeeeee", xpd=T)

i_stop = 40
y = 33.4
xi= 0.5
yi = 6.8/length(rows)
palette = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(101))
for (row in rows) {
  for (i in 1:80) {
    x = i 
    if (i > i_stop) { x = x + 2 }
    col = round((tcga_classifier[row, i] - min(tcga_classifier[row, ])) / (max(tcga_classifier[row, ]) - min(tcga_classifier[row, ])) * 100) + 1
    rect(x-xi, y-yi/2, x+xi, y+yi/2, border=NA, col=palette[col],xpd=T)
  }
  y = y - yi
}

# Defaults
white = "#EEEEEE"
border = NA
xi = 0.4

# GEP
y = 26
for (i in 1:80) {
  col = white
  x = i 
  if (i > i_stop) { x = x + 2 }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[[tcga_data$GEP[i]]])
}

# Gaq mutations
y = 25
for (gaq in c("GNAQ", "GNA11", "CYSLTR2", "PLCB4")) {
  y = y-1
  for (i in 1:80) {
    x = i 
    if (i > i_stop) { x = x + 2 }
    val = tcga_data$`Gαq mutation`[i]
    col = white
    if (!is.na(val) & stringr::str_detect(val, gaq)) {
      col = "#777777"
    }
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  }
}

# EIF1AX
y = 17
for (i in 1:80) {
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`EIF1AX mutation`[i]
  col = white
  if (!is.na(val)) {
    col = "#777777"
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}

# SF3B1
y = 18
for (i in 1:80) {
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`SF3B1 mutation`[i]
  col = white
  if (!is.na(val)) {
    col = "#777777"
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}

# BAP1
y = 19
for (i in 1:80) {
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`BAP1 mutation`[i]
  col = white
  if (!is.na(val)) {
    col = "#777777"
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}

# Chr. 3
y = 15
for (i in 1:80) {
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`Chr. 3 (DNA)`[i]
  col = white
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  if (!is.na(val)) {
    col = white
    if (val == "trisomy") { 
      col = "#AAAAAA"
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    } 
    else {
      if (val == "monosomy") { 
        col = "#777777"
        rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
      }
      else if (val == "isodisomy") { 
        col = "#777777"
        rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
      }
      else if (val == "monosomy (arm-level variation)") { 
        col = "#F44336"
        rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col="#777777")
      } 
      
    }
  }
}

# Chr. 8q
y = 14
for (i in 1:80) {
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`Chr. 8q (DNA)`[i]
  col = white
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
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
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    
    if (val %in% c("gain/amplification", "gain/amplification (biallelic)")) {
      col = "#777777"
      rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
      col = "#AAAAAA"
      polygon(c(x-xi, x-xi, x+xi), c(y-0.4, y+0.4, y+0.4), border=border, col=col)
    }
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
y = 7+4
x =xn+1
text(x-2.5, y, labels = "GEP classifier", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#AAAAAA", xpd=T)
text(x, y, labels = "class I", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#777777", xpd=T)
text(x, y, labels = "class II", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

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

y=y-3
x=xn+1
text(x-2.5, y, labels = "Mutations", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white, xpd=T)
text(x, y, labels = "wild-type", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn+18.5

rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#777777", xpd=T)
text(x, y, labels = "mutant", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


y = yn
y = 7+4
xn= 82-39
x=xn
text(x-2.5, y, labels = "Chromosome 3", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#777777", xpd=T)
text(x, y, labels = "loss", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn+15
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white, xpd=T)
rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col="#777777", xpd=T)

text(x, y, labels = "partial loss", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn

y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white, xpd=T)
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#AAAAAA", xpd=T)
text(x, y, labels = "gain", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


y=y-2
text(x-2.5, y, labels = "Chromosome 8q", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#AAAAAA", xpd=T)
text(x, y, labels = "gain", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn+15
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#777777", xpd=T)
text(x, y, labels = "amplification", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

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

# Plot and save dendrogram
file = "gene-expression-profiling/figure-1C-dendrogram.png"
png(file, res=600, width=5000, height=1500)
dendrogram = pheatmap(tcga_classifier, 
                      scale="row", 
                      custering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                      cutree_cols = 2, 
                      cutree_rows = 2,
                      treeheight_col=35)
dev.off()
system(paste("open", file))

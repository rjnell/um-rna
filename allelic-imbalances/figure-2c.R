###
##  Figure 2C
###

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

# LUMC-data
tcga_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=2)


tcga_data = tcga_data[order(tcga_data$Order,decreasing = F),]
#tcga_data = tcga_data[order(tcga_data$GEP),]

#tcga_data = tcga_data[match(lumc_classifier_ph$tree_col$labels[lumc_classifier_ph$tree_col$order], tcga_data$Tumor_ID),]

colors = list()
colors[["female"]] = RColorBrewer::brewer.pal(12, "Paired")[5]
colors[["Class I"]] = "#79BFF7"
colors[["Class II"]] = RColorBrewer::brewer.pal(12, "Paired")[5]
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

file = "allelic-imbalances/figure-2c.png"
png(file, res=600, width=5000, height=7500)
par(mar=c(5,8,5,3))
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

y = 25
white = "#DDDDDD"
border = NA
xi = 0.4

y = 24
for (i in 1:80) {
  col = white
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`Chr. 3 (DNA)`[i]
  if (val != "disomy") {
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    col = '#999999'
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  }
  else {
    col = white
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  }
}

y = 21
for (i in 1:80) {
  col = white
  x = i 
  if (i > i_stop) { x = x + 2 }
  val = tcga_data$`Chr. 8q (RNA)`[i]
  if (val == "LOH") {
    col = '#999999'
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
y = 23
for (i in 1:80) {
  col = white
  x = i 
  if (i > i_stop) { x = x + 1 }
  val = tcga_data$`Chr. 3 (DNA)`[i]
  val_detect = tcga_data$`Chr. 3 (DNA)`[i]
  col = "#DDDDDD"  
  if (val_detect == "monosomy (arm-level variation)") { 
    col = "#DDDDDD"
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    col = "#2196F3"
    rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col=col)
  } else {
    if (val_detect == "disomy") { 
      
      
    }
    else if (val_detect == "isodisomy") { 
      col = "#0B72C4"
    }
    else if (val_detect == "trisomy") { 
      col = "#F44336"
      
    } 
    else { 
      col = "#2196F3"
      
    }
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  }
}
y = 20
for (i in 1:80) {
  col = white
  x = i 
  if (i > i_stop) { x = x + 1 }
  val = tcga_data$`Chr. 8q (DNA)`[i]
  col = "#DDDDDD"  
  if (val == "gain") { 
    col = "#F44336"
  } 
  else if (val == "amplification (biallelic)") { 
    col = "darkorange"
  } 
  else if (val == "amplification") { 
    col = "brown"#"#D4190C"
  } 
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
  
  if (val == "gain/amplification") {
    col = "brown"#"#777777"
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    col = "#F44336"
    polygon(c(x-xi, x-xi, x+xi), c(y-0.4, y+0.4, y+0.4), border=border, col=col)
  }
  if (val == "gain/amplification (biallelic)") {
    col = "darkorange"#"#777777"
    rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
    col = "#F44336"
    polygon(c(x-xi, x-xi, x+xi), c(y-0.4, y+0.4, y+0.4), border=border, col=col)
  }
}



pos = 1-0.7-2
axis(side = 2, at = 24:23, labels = c("RNA", "DNA"), las=2, 
     cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
text(-13, 23.5, labels = "Chromosome 3", cex=1.1, col="#333333", font=1, pos=2, xpd=T)
axis(side = 2, at = 21:20, labels = c("RNA", "DNA"), las=2, 
     cex.axis=1.1, col="#333333", pos = pos, col.ticks = "#b1b1b1", lwd.ticks = 1.4)
text(-13, 20.5, labels = "Chromosome 8q", cex=1.1, col="#333333", font=1, pos=2, xpd=T)
segments(pos, 24+0.4, pos, 23-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 21+0.4, pos, 20-0.4, lwd=1.4, col="#b1b1b1")

xn = -22.5
y = 9
x =xn

yn = 16+1.5
y = yn
xn = 1
x = xn
text(x-2.5, y, labels = "RNA", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["DNA + RNA"]])
text(x, y-0.05, labels = "LOH", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y=y-1
x=xn
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white)
text(x, y-0.05, labels = "no LOH", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


yn = 16+1.5
y = yn
xn = 30
x = xn
text(x-2.5, y, labels = "DNA", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white,xpd=T)
text(x, y-0.05, labels = "disomy", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=x+23
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white, xpd=T)
rect(x-xi, y-0.2, x+xi, y+0.2, border=border, col="#2196F3", xpd=T)
text(x, y-0.05, labels = "partial loss", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn
y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#2196F3", xpd=T)
text(x, y-0.05, labels = "loss", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=x+23
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#0B72C4", xpd=T)
text(x, y-0.05, labels = "isodisomy", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn
y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#F44336", xpd=T)
text(x, y-0.05, labels = "gain", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=x+23
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="brown", xpd=T)
text(x, y-0.05, labels = "amplification", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn
y=y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="darkorange", xpd=T)
text(x, y-0.05, labels = "amplification (biallelic)", cex=1.1, col="#333333", font=1, pos=4, xpd=T)

dev.off()
system(paste("open", file))





###
##  Figure 1 - Landscape of transcriptional, molecular and clinical characteristics
###

# Close current images
system("taskkill /F /IM PhotosApp.exe /T")
graphics.off()

# Load libraries
library(readxl)
library(stringr)

tcga_data = read_xlsx("data/Supplementary Tables.xlsx", sheet=2)

tcga_data = tcga_data[order(tcga_data$`Gaq detectability`),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "PLCB4")),which(!str_detect(tcga_data$`Gαq mutation`, "PLCB4")),which(is.na(tcga_data$`Gαq mutation`))),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "CYSLTR2")),which(!str_detect(tcga_data$`Gαq mutation`, "CYSLTR2")),which(is.na(tcga_data$`Gαq mutation`))),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "GNA11")),which(!str_detect(tcga_data$`Gαq mutation`, "GNA11")),which(is.na(tcga_data$`Gαq mutation`))),]
tcga_data = tcga_data[c(which(str_detect(tcga_data$`Gαq mutation`, "GNAQ")),which(!str_detect(tcga_data$`Gαq mutation`, "GNAQ")),which(is.na(tcga_data$`Gαq mutation`))),]

tcga_data = tcga_data[order(tcga_data$`BAP1 detectability`),]
tcga_data = tcga_data[order(tcga_data$`SF3B1 detectability`),]
tcga_data = tcga_data[order(tcga_data$`EIF1AX detectability`),]

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

file = "mutations/figure-3A.png"
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


y=25
for (gaq in c("GNAQ", "GNA11", "CYSLTR2", "PLCB4")) {
  # SF3B1
  y = y-1
  for (i in 1:80) {
    x = i 
    if (i > i_stop) { x = x + 2 }
    val = tcga_data$`Gαq mutation`[i]
    col = white
    if (!is.na(val) & stringr::str_detect(val, gaq)) {
      det = tcga_data$`Gαq detectability`[i]
      if (det == "DNA") {
        col = colors[["DNA"]]
      }
      else if (det == "RNA") {
        col = colors[["RNA"]]
      }
      else {
        col = "#777777"  
      }
    }
    
    if (tcga_data$`Patient ID`[i] == "TCGA-YZ-A985") {
      if (gaq == "GNAQ") {
        col = colors[["DNA"]]
      }
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
    det = tcga_data$`EIF1AX detectability`[i]
    if (det == "DNA") {
      col = colors[["DNA"]]
    }
    else if (det == "RNA") {
      col = colors[["RNA"]]
    }
    else {
      col = "#777777"  
    }
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
    det = tcga_data$`BAP1 detectability`[i]
    if (det == "DNA") {
      col = colors[["DNA"]]
    }
    else if (det == "RNA") {
      col = colors[["RNA"]]
    }
    else {
      col = "#777777"  
    }
  }
  rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=col)
}
pos = 1-0.7-2

limits=c(33,26)
#segments(pos, limits[1]+0.4, pos, limits[2]-0.4, lwd=1.4, col="#b1b1b1", xpd=T)

axis(side = 2, at = 24:21, labels = c("GNAQ", "GNA11", "CYSLTR2", "PLCB4"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)



axis(side = 2, at = 17:19, labels = c("EIF1AX","SF3B1","BAP1"), las=2, cex.axis=1.1, col="#333333", pos = pos, font=3, col.ticks = "#b1b1b1", lwd.ticks = 1.4)

#pos = c(pos, s)
#segments(pos, 26+1.4, pos, 26-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 21+3.4, pos, 21-0.4, lwd=1.4, col="#b1b1b1")
segments(pos, 17+2.4, pos, 17-0.4, lwd=1.4, col="#b1b1b1")

y = yn
y = 24
xn= 90
x=xn
text(x-2.5, y, labels = "Detectability", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col="#777777", xpd=T)
text(x, y, labels = "DNA + RNA", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y =y -1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["DNA"]], xpd=T)
text(x, y, labels = "DNA only", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
x=xn

y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=colors[["RNA"]], xpd=T)
text(x, y, labels = "RNA only", cex=1.1, col="#333333", font=1, pos=4, xpd=T)
y = y-1
rect(x-xi, y-0.4, x+xi, y+0.4, border=border, col=white, xpd=T)
text(x, y, labels = "not detected", cex=1.1, col="#333333", font=1, pos=4, xpd=T)


dev.off()
system(paste("open", file))


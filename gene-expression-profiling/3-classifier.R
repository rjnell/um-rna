###  
### STEP 3 - GET CLASSIFIER FROM TCGA
###

# Set seed
set.seed(83539)

# Load normalised TCGA data
tcga_normalised = as.matrix(read.csv("gene-expression-profiling/tcga-normalised.tsv",sep="\t", stringsAsFactors = F, check.names=F, row.names = 1))

# Load TCGA molecular data
molecular_data = readxl::read_excel("data/Supplementary Tables.xlsx", sheet=2)
molecular_data = molecular_data[match(colnames(tcga_normalised),molecular_data$`Patient ID`),]

# Obtain twoclass mRNA categorization
tcga = tcga_normalised[,order(colnames(tcga_normalised))]
molecular_data = molecular_data[order(molecular_data$`Patient ID`),]
twoclass = molecular_data$GEP

mol = cbind(molecular_data, twoclass)
rownames(mol) = molecular_data$`Patient ID`
colnames(tcga)==rownames(mol)

# Load ClaNC
source("bin/clanc.R")

# Train and build clanc for two-class signature, with 100 genes per class
no_of_classes = 2
tumor_classes = mol$twoclass 
trainClaNC_data = trainClanc(tcga, tumor_classes, rownames(tcga))
n=100
buildClaNC_data = buildClanc(tcga, tumor_classes, 1:no_of_classes, trainClaNC_data, active = n)

# Annotate and save classifier signature
classifier = matrix(data="",nrow=n*2,ncol=2)
colnames(classifier) = c("ENSG", "class")
classifier[,1] = buildClaNC_data$geneNames
classifier[which(buildClaNC_data$cntrds[,1]>buildClaNC_data$cntrds[,2]),2] = "Class I higher"
classifier[which(buildClaNC_data$cntrds[,2]>buildClaNC_data$cntrds[,1]),2] = "Class II higher"
ensg_annotation = read.csv("data/ensg-annotation.tsv", sep="\t", check.names = F, stringsAsFactors = F)
classifier_annotated = cbind(classifier, ensg_annotation[match(classifier[,1], ensg_annotation$ensembl_gene_id),])
classifier_annotated[,1] == classifier_annotated[,3] 
write.table(classifier_annotated[,c(1:2,4:7)], "gene-expression-profiling/classifier.tsv", sep="\t")

View(classifier_annotated[,c(1:2,4:7)])
# Save image of classifier in TCGA data
tcga_clanc = tcga_normalised[as.character(classifier_annotated[,1]),]
tcga_clanc = tcga_clanc[which(matrixStats::rowSds(tcga_clanc)>0.5,),]
ann_tcga = data.frame(mol[,c(3,4,5)])
png("gene-expression-profiling/classifier-tcga.png",res=600,width=5000,height=5000)
pt = pheatmap::pheatmap(tcga_clanc, 
                        annotation_col = ann_tcga, 
                        scale="row", 
                        clustering_callback = function(hc, ...) { dendsort::dendsort(hc, type = "min") }, 
                        cutree_cols = 2, 
                        cutree_rows = 2, 
                        color = RColorBrewer::brewer.pal(9,"RdBu"), 
                        show_rownames = T, 
                        fontsize_col = 5,
                        fontsize_row = 2,
                        breaks = seq(from=-7, to=7, length.out = 10))
dev.off()




pie = function (f1,f2,limit, x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", 
        "lightcyan", "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab) & dx[i]>limit) {
      lines(c(1, f1) * P$x, c(1, f1) * P$y, col = "#b1b1b1", lwd=1.4)
      text(f2 * P$x, f2 * P$y, labels[i], xpd = TRUE, col="#333333",
           , ...) #adj = ifelse(P$x < 0, 1, 0)
    }
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    
  }
  title(main = main, ...)
  invisible(NULL)
}


chromosomes = c(1:22, "X", "Y")
res = matrix(0, ncol=24, nrow=2)
colnames(res) = chromosomes

t = table(classifier_annotated$chromosome_name[which(classifier_annotated$class == "Class I higher")])
for (i in 1:length(t)) {
  res[1,names(t)[i]] = t[i]
}

t = table(classifier_annotated$chromosome_name[which(classifier_annotated$class == "Class II higher")])
for (i in 1:length(t)) {
  res[2,names(t)[i]] = t[i]
}


colors = colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))(24)
names(colors) = chromosomes

png("gene-expression-profiling/figure-1B1.png", res=600, width=2500, height=2500)
cl1 = res[1,]
sel = which(cl1 != 0)
pie(cl1[sel], col=colors[sel],init.angle = 90,clockwise = T, border=F, limit = 4/165,f1=1.15, f2=1.3)
dev.off()

png("gene-expression-profiling/figure-1B2.png", res=600, width=2500, height=2500)
cl2 = res[2,]
sel = which(cl2 != 0)
pie(cl2[sel], col=colors[sel],init.angle = 90,clockwise = T, border=F, radius = sqrt(35/165 / pi) * 0.8, limit = 2/35,f1=1.5,f2=2)
dev.off()







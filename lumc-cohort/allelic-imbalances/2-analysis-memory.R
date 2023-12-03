
sample_id = "GS_103114-001-032"

# Set output paths
output_path = "/mnt/d/SURFDRIVE/ORGANIZED/R/um-rna/lumc-cohort/allelic-imbalances/raw/"
output_path_win = "D:/SURFDRIVE/ORGANIZED/R/um-rna/lumc-cohort/allelic-imbalances/raw/"

# Create sample_id-specific output dir
output_path = paste0(output_path, sample_id, "/")
output_path_win = paste0(output_path_win, sample_id, "/")

# Annotate variants
annotation_file = "/mnt/f/dbsnp151.vcf.gz"
annotation_cmd = paste0("wsl java -jar /mnt/d/SURFDRIVE/ORGANIZED/R/um-rna/bin/snpEff/snpSift.jar annotate ", annotation_file, " ", 
                        paste0(output_path, sample_id, ".recode.chr.vcf.gz"), " > ", 
                        paste0(output_path, sample_id, "-annotated.vcf"))
system(annotation_cmd)


# Read vcf and prepare table
vcf_file = vcfR::read.vcfR(paste0(output_path_win, sample_id, "-annotated.vcf"))      
cols = c("GT","ABQ","AD","ADF","ADR","DP","FREQ","GQ","PVAL","RBQ","RD","RDF","RDR","SDP","COMMON","HLA")
vcf_table = cbind(rep(sample_id, nrow(vcf_file@fix)), rep("gene_target", nrow(vcf_file@fix)), vcf_file@fix, matrix(data=NA, nrow=nrow(vcf_file@fix), ncol=length(cols)))
colnames(vcf_table) = c("sample", "RNA-seq filtered", colnames(vcf_file@fix), cols)
for (index in row(vcf_file@fix)[,1]) {
  formats = strsplit(vcf_file@gt[index, 1], split=":")[[1]]  
  values = strsplit(vcf_file@gt[index, 2], split=":")[[1]]  
  vcf_table[index,formats] = values
  vcf_table[index,"COMMON"] = grepl("COMMON=1", vcf_file@fix[index, "INFO"], fixed=T)
  vcf_table[index,"HLA"] = grepl("HLA-", vcf_file@fix[index, "INFO"], fixed=T)
}
head(vcf_table)

s = 1
vcf_selected = vcf_table[which(!is.na(vcf_table[,"AD"]) & vcf_table[,"AD"]!="."),]
vcf_selected = vcf_selected[which(as.numeric(vcf_selected[,"AD"]) > 2*s),]
vcf_selected = vcf_selected[which(as.numeric(vcf_selected[,"RD"]) > 2*s),]
vcf_selected = vcf_selected[which(as.numeric(vcf_selected[,"AD"])+as.numeric(vcf_selected[,"RD"]) > 40*s),]
vcf_selected = vcf_selected[which(vcf_selected[,"COMMON"] == "TRUE"),]
vcf_selected = vcf_selected[which(vcf_selected[,"HLA"] == "FALSE"),]
vcf_selected = vcf_selected[order(as.numeric(vcf_selected[,"POS"])),]
vcf_selected[which(vcf_selected[,"CHROM"] == "X"),"CHROM"]=23
vcf_selected = vcf_selected[order(as.numeric(vcf_selected[,"CHROM"])),]
vcf_selected[which(vcf_selected[,"CHROM"] == 23),"CHROM"]="X"

# VAFS
vafs = as.numeric(sub("%", "", vcf_selected[,"FREQ"]))

# Our chromosomes
chromosome_list = c(1:22, "X")
chromosome_data = table(vcf_selected[,"CHROM"])
chromosome_data = chromosome_data[match(chromosome_list, names(chromosome_data))]

# Centromeres
centromere_data = read.csv("bin/centromeres.txt", stringsAsFactors = F, header=F, sep="\t")
centromeres = NULL
mins = NULL
for (chr in chromosome_list) {
  min = max(centromere_data$V3[which(centromere_data$V2 == paste0("chr",chr))])
  centromeres = c(centromeres, min)
  sel = which(vcf_selected[,"CHROM"] == chr & as.numeric(vcf_selected[,"POS"]) < min)
  if (length(sel) > 0) {
    mins = c(mins, max(sel))
  }
}
names(centromeres) = c(1:22,"X")

# Summarise results
baf = NULL
for (i in 1:nrow(vcf_selected)) {
  b = vcf_selected[i,c("CHROM","POS","FREQ", "AD", "RD", "DP")]
  b[3] = sub("%", "", b[3])
  if (as.numeric(b[2]) < centromeres[b[1]]) {
    b = c(b,paste0(b[1], "p"))
  }
  else {
    b = c(b,paste0(b[1], "q"))
  }
  baf = rbind(baf, b)
}
colnames(baf) = c("CHR", "POS", "BAF", "AD", "RD", "DP", "ARM")
rownames(baf) = NULL
write.table(baf, paste0(output_path_win, sample_id,"_baf.tsv"), sep="\t")

# Filename of image
result_file = paste0(output_path_win, sample_id,".png")

# Generate image
graphics.off()
png(result_file, res=600, width=5500, height=1800)  

par(mar=c(5,5,5,5))
ylim = c(0,1)
xlim = c(0,20)
plot(type="n",
     x = xlim,
     y = ylim,
     axes = F,
     ylab="",
     xlab="",
     main = sample_id,
     bty = "l", 
     xaxs = "i", 
     yaxs = "i")

yat = c(0,1/3,1/2,2/3,1)#seq(from=ylim[1],to=ylim[2],by=0.1)
segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],xpd=T,col="#B1B1B1",lwd=1.4)

# Y axis
axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.053,lwd=1.4)
axis(side = 2,at=yat,las=2,labels=c("0%","33%","50%","67%","100%"), lwd=0, col.axis="#333333",line=-0.23)
mtext(text = "B-allele fraction", side = 2, line=3.8,col="#333333")  

scalefactor = 1/length(vafs)*20

for (min in mins) {
  segments(min*scalefactor,0,min*scalefactor,1,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
}

pos = 0.5
for (chr in 1:length(chromosome_data)) {
  chromosome = names(chromosome_data)[chr]
  pos_new = pos + chromosome_data[chr]  
  segments(pos_new*scalefactor,0,pos_new*scalefactor,1,col = "#b1b1b1",lwd=1.4,xpd=T)
  #mtext(at = (pos + pos_new)/2, line = 0, side = 1, chromosome, cex=0.5)
  
  text((pos + pos_new)/2*scalefactor, -0.1, srt=60, adj=1, col="#333333", labels=paste0(chr), xpd=T)  
  pos = pos_new
}

points((1:length(vafs))*scalefactor,vafs/100,pch = 16,
       cex = 0.23,col=rgb(0,0,0,0.2))


dev.off()
system(paste0("open ", result_file))
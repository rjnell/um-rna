###
##  Script to generate figures from RNA and DNA data
###

# Smoothing parameter range [x-n:x+n]
n = 99

# Number of 
g_size = 10
s_limit = 25

# Function to perform segmentation
segmentation = function(bafs, n=49, limit=40, g_size=20) {
  
  #bafs = baf$BAF[which(baf$CHR == 17)]
  #bafs = c(bafs)
  #plot(1:length(bafs),bafs)
  
  # Mirror BAFs
  bafs[which(bafs>50)] = 100-bafs[which(bafs>50)]
  
  # Initiate BAF groups
  groups = seq(from=0, to=50, by=50/g_size) #seq(from=-1.25, to=51.25, by=2.5)
  
  # Initiate vector of medians
  m = NULL
  
  # Initiate vector of segment ranges
  s = 1
  
  # Iterate through values
  for (x in 1:length(bafs)) {
    
    # Determine range based on smoothing parameter n
    range = max(1,x-n):min(length(bafs),x+n)
    
    # Summarize values according to groups
    gs = NULL
    for (g in 1:(length(groups)-1)) {
      gs[g] = length(which(bafs[range] >= groups[g] & bafs[range] < groups[g+1]))
    }
    mm = groups[which(gs==sort(gs, TRUE)[1])[1]]+25/g_size
    mm = median(bafs[range])
    m = c(m, mm)
    
    # Check if through 40-value ==> segment range
    if (x>1) {
      if (m[x] < limit & m[x-1] > limit) {
        s = c(s, min(length(bafs), round(x+length(range)/4)))
      }
      if (m[x] > limit & m[x-1] < limit) {
        s = c(s, max(1, round(x-length(range)/4)))
      }
    }
    #points(x,m[x],pch=16,col="blue")
    
  }
  
  # Finalize vector of segment ranges
  s = c(s, length(bafs))
  s = unique(s)
  
  #segments(s,0,s,100)
  
  # Initialize output segments
  ss = NULL
  
  # Iterate through segment ranges
  for(i in 1:(length(s)-1)) {
    
    # Select values in segment
    v = bafs[s[i]:s[i+1]]
    
    # Summarize values according to groups
    gs = NULL
    for (g in 1:(length(groups)-1)) {
      gs[g] = length(which(v >= groups[g] & v < groups[g+1]))
    }
    
    # Determine largest group and correct for noise
    val = groups[which(gs==sort(gs, TRUE)[1])[1]]+25/g_size
    #val = median(v)
    if (val > limit) { val = 50 }
    else { val = median(v) }
    
    # Save results
    ss = rbind(ss, c(s[i], s[i+1], val, 100-val))
    
  }
  
  # Return results
  return(ss)
}

# Read samplesheet data
samplesheet_original = read.table("data/tcga-samplesheet.tsv", sep="\t")[2:81,]
samplesheet = cbind(substr(samplesheet_original$V2,0,36), samplesheet_original$V6)

# Read raw files
files = dir("allelic-imbalances/raw/")
files = files[which(substr(files,0,2) == "s_")]
files = substr(files,3,nchar(files))
files = gsub("_","-",files)

# Read (allele-specific) copy number data from DNA
cn_data = data.table::fread("data/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg", sep="\t")
allele_specific_cn_data = data.table::fread("../copy-number-alterations/data/TCGA_mastercalls.abs_segtabs.fixed.txt", sep="\t")

# Iterate through samples
for (s_id in which(samplesheet[,1] %in% files)[c(67,52)]) {
  
  # Load sample/ID information
  sample = gsub("-","_",samplesheet[s_id,1])
  TCGA_id = samplesheet[s_id,2]
  
  # Load BAF information
  baf = read.table(paste0("allelic-imbalances/raw/s_", sample, "/s_", sample, "_baf.tsv"))
  
  # Exclude X chromosome
  baf = baf[which(baf$CHR != "X"),]
  baf = cbind(baf, rep(NA, nrow(baf)), rep(NA, nrow(baf)))
  colnames(baf)[8:9] = c("CN", "LOH")
  
  # Load sample-specific copy number data
  uvm_cn_data = data.frame(cn_data[which(substr(cn_data$Sample,0,15) 
                                         %in% paste0(TCGA_id,"-01")),])
  uvm_loh_data = data.frame(allele_specific_cn_data[which(allele_specific_cn_data$Sample
                                                          %in% paste0(TCGA_id,"-01")),])
  
  # Iterate through BAFs and link copy number data
  for (j in 1:nrow(baf)) {
    chr = baf$CHR[j]
    if (chr == "X") { chr = 23 }
    w = which(uvm_cn_data$Chromosome == as.numeric(chr) &
                uvm_cn_data$Start < as.numeric(baf$POS[j]) &
                uvm_cn_data$End > as.numeric(baf$POS[j]))
    w_loh = which(uvm_loh_data$Chromosome == as.numeric(chr) &
                    uvm_loh_data$Start < as.numeric(baf$POS[j]) &
                    uvm_loh_data$End > as.numeric(baf$POS[j]))
    if (length(w) == 1 & length(w_loh) == 1) {
      baf$CN[j] = uvm_cn_data$Segment_Mean[w]
      baf$LOH[j] = uvm_loh_data$LOH[w_loh]
    }
  }
  
  # Only keep BAFs with complete copy number data
  baf = baf[which(!is.na(baf$CN) & !is.na(baf$LOH)),]
  
  # Figure with manual segmentation
  result_file = paste0("allelic-imbalances/res/manual-segmentation/", TCGA_id, ".png")
  png(result_file, res=600, width=7000, height=2850)  
  {
    par(mar=c(5,7,5,3))
    ylim = c(-1.5,1.2)
    xlim = c(0,20)
    plot(type="n",
         x = xlim,
         y = ylim,
         axes = F,
         ylab="",
         xlab="",
         bty = "l", 
         xaxs = "i", 
         yaxs = "i")
    
    ym = 0.2
    chrmove=-0.25
    yat = c(0,1/3,1/2,2/3,1)+ym
    segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
    segments(xlim[1],min(yat),xlim[1],max(yat),xpd=T,col="#B1B1B1",lwd=1.4)
    
    # Y axis
    axis(side = 2,at=yat,las=2,labels=c("0%","33%","50%","67%","100%"),col.ticks = "#b1b1b1", 
         col = "#b1b1b1", col.axis="#333333", lwd=1.4,cex=1.1,cex.axis=1.1)
    mtext(text = "BAF", side = 2, line=4.5, at=0.5+ym, col="#333333",cex=1.1,cex.axis=1.1)  
    mtext(text = "(RNA)", side = 2, line=3.5, at=0.5+ym, col="#333333",cex=1.1,cex.axis=1.1)
    
    
    # Comparison:
    ymax = max(c(1,ceiling(quantile(2^baf$CN*2,0.99))-2))
    yat = c(1,2:(2+ymax))/6-(5+ymax)/6
    segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
    #segments(xlim[1],0,xlim[1],ylim[2],xpd=T,col="#B1B1B1",lwd=1.4)
    axis(side = 2,at=yat,las=2,labels=c("\u22121","0",paste0("+",1:ymax)),
         col.ticks = "#b1b1b1", col = "#b1b1b1", col.axis="#333333", lwd=1.4,cex=1.1,cex.axis=1.1)
    axis(side = 2,at=min(yat)-.15-1/12,las=2,labels=c("LOH"),col.ticks = "#b1b1b1", 
         col = "#b1b1b1", col.axis="#333333", lwd=1.4,cex=1.1,cex.axis=1.1)
    axis(side = 2,at=chrmove+0.1,las=2,labels=c("Chromosome"),col.ticks = "#b1b1b1", 
         col = "#b1b1b1", col.axis="#333333", lwd=1.4,cex=1.1,cex.axis=1.1)
    
    mtext(text = "Copy", side = 2, line=5.5, at=mean(c(min(yat)-.15-1/6, -0.5)), col="#333333",cex=1.1,cex.axis=1.1) 
    mtext(text = "number", side = 2, line=4.5, at=mean(c(min(yat)-.15-1/6, -0.5)), col="#333333",cex=1.1,cex.axis=1.1) 
    mtext(text = "(DNA)", side = 2, line=3.5, at=mean(c(min(yat)-.15-1/6, -0.5)), col="#333333",cex=1.1,cex.axis=1.1)
    
    # Define scalefactor (depending on total number of BAFs)
    scalefactor = 1/nrow(baf)*20
    
    # Plot copy number value & LOH status
    for (b in 1:nrow(baf)) {
      
      # Copy number value
      sm = baf$CN[b]
      cnv = 2^sm*2
      
      if (sm < -0.3) {
        if (cnv > 1.6) {
          col = "#A6D5FA"
        } else {
          col = "#2196F3" ##0B72C4"
        }
        
      }
      else if (sm > 0.3) {
        
        if (cnv > 3) {
          col = "#D4190C"
        }
        else if (cnv < 2.4) {
          col = "#FBB4AF"
        }
        else {
          col = "#F44336"
        }
        
      }
      else {
        col = "#333333"
      }
      #segments(b*scalefactor, -0.5, b*scalefactor, -0.6, col=col)
      #points(b*scalefactor, cnv/6-7/6, pch=16, cex=0.23, xpd=T, col= col)
      if (col != "#333333") {
        segments(b*scalefactor, 2/6-(5+ymax)/6, b*scalefactor, cnv/6-(5+ymax)/6, xpd=T, col= col)
      } 
      
      if (baf$LOH[b] == 1) {
        col = "#333333"
      }
      else {
        col = "#ffffff"
      }
      segments(b*scalefactor, min(yat)-0.15, b*scalefactor, min(yat)-.15-1/6, col=col, xpd=T)
    }
    
    #segments(0,c(min(yat)-0.15,min(yat)-0.275),20,lwd=1.4,col="#eeeeee",xpd=T)
    # Dotted lines between each chromosomal arm
    for (arm in unique(baf$ARM)) {
      min = min(which(baf$ARM == arm))
      segments(min*scalefactor,0+ym,min*scalefactor,1+ym,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
      segments(min*scalefactor,min(yat),min*scalefactor,-0.5,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
      segments(min*scalefactor,min(yat)-0.15,min*scalefactor,min(yat)-0.15-1/6,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
    }
    
    # Solid lines between each chromosome + labels
    chromosome_data = table(baf$CHR)
    chromosome_data = chromosome_data[match(1:22, names(chromosome_data))]
    pos = 0.5
    for (chr in 1:length(chromosome_data)) {
      chromosome = names(chromosome_data)[chr]
      pos_new = pos + chromosome_data[chr]  
      
      segments(pos_new*scalefactor,0+ym,pos_new*scalefactor,1+ym,col = "#b1b1b1",lwd=1.4,xpd=T)
      segments(pos_new*scalefactor,min(yat),pos_new*scalefactor,-.5,col = "#b1b1b1",lwd=1.4,xpd=T)
      segments(pos_new*scalefactor,min(yat)-0.15,pos_new*scalefactor,min(yat)-0.15-1/6,col = "#b1b1b1",lwd=1.4,xpd=T)
      if (chr%%2 == 0) {
        #rect(pos*scalefactor, chrmove,pos_new*scalefactor,chrmove+0.2,col="#eeeeee",border=NA,xpd=T) 
      }
      segments(pos_new*scalefactor,chrmove,pos_new*scalefactor,chrmove+0.2,col = "#b1b1b1",lwd=1.4,xpd=T)
      segments(pos*scalefactor,chrmove,pos*scalefactor,chrmove+0.2,col = "#b1b1b1",lwd=1.4,xpd=T)
      
      if ((pos_new - pos)*scalefactor > 0.45 | chr<10) {
        text((pos + pos_new)/2*scalefactor, chrmove+0.11, col="#333333", labels=paste0(chr), xpd=T, cex=1.1) 
      }
      else {
        if (chr != 22) {
          if (chromosome_data[chr-1] > sum(chromosome_data)/20) {
            #if (chromosome_data[chr+1] > sum(chromosome_data)/20) {
            #text((pos + pos_new)/2*scalefactor, 1.11, col="#333333", labels=paste0(chr), xpd=T, cex=1.1) 
            #}
          }
        }
      }
      #text((pos + pos_new)/2*scalefactor, -0.26, col="#333333", labels=paste0(chr), xpd=T, cex=1.1) 
      pos = pos_new
    }
    
    # Plot BAFs
    points((1:nrow(baf))*scalefactor,baf$BAF/100+ym,pch = 16,
           cex = 0.3,col=rgb(0,0,0,0.2),xpd=T)
    
    # Overlay segments
    segments(0,min(yat)-0.15,0,min(yat)-0.15-1/6,lwd=1.4,col="#b1b1b1",xpd=T)
    segments(0,2/6-(5+ymax)/6,20,col="#b1b1b1",lwd=1.4)
    #segments(0,chrmove,20,col="#eeeeee",lwd=1.4)
    #segments(0,chrmove+0.2,20,col="#eeeeee",lwd=1.4)
    segments(0,chrmove+0.2,0,chrmove,col="#b1b1b1",lwd=1.4)
    
    # Title
    y=1.5
    text(-0.2, y, col="#333333", labels=TCGA_id, xpd=T, cex=1.1, font=2, pos=4) 
    
    # Labels
    x = 7.1
    text(x+0.3,y, cex=1.1,col="#333333", labels="BAF segmentation",xpd=T, pos=4)
    xx = c(0.08,0.18,0.00,0.03,0.23,0.13,0.09,0.23,0.12,0.12,0.26,0.09,0.18,0.10,0.14,0.20,0.14,0.10,0.03,0.26,0.21,0.26,0.16,0.13,0.13,0.17,0.21,0.07,0.12,0.19)
    yy = c(1.51,1.51,1.52,1.52,1.51,1.46,1.48,1.48,1.52,1.51,1.47,1.48,1.50,1.51,1.54,1.50,1.49,1.50,1.47,1.49,1.51,1.52,1.47,1.47,1.51,1.47,1.51,1.55,1.51,1.51)
    points(x+xx+0.01,yy,pch = 16,
           cex = 0.3,col=rgb(0,0,0,0.2),xpd=T)
    segments(x,y,x+0.3, lwd=1.4, col = "#333333", xpd=T)
    
    x = 11.8
    rect(x,y-1/12,x+.1,y+1/12, col="#F44336", border=NA, xpd=T)
    text(x+0.1,y, cex=1.1,col="#333333", labels="gain",xpd=T, pos=4)
    
    x = 13.8
    rect(x,y-1/12,x+.1,y+1/12, col="#D4190C", border=NA, xpd=T)
    text(x+0.1,y, cex=1.1,col="#333333", labels="amplification",xpd=T, pos=4)
    
    x = 17.1
    rect(x,y-1/12,x+.1,y+1/12, col="#2196F3", border=NA, xpd=T)
    text(x+0.1,y, cex=1.1,col="#333333", labels="loss",xpd=T, pos=4)
    
    x = 18.95
    rect(x,y-1/12,x+.1,y+1/12, col="#333333", border=NA, xpd=T)
    text(x+0.1,y, cex=1.1,col="#333333", labels="LOH",xpd=T, pos=4)
  }
  for (arm in unique(baf$ARM)) {
    
    # Load arm data
    arm_data = baf[which(baf$ARM == arm),]
    
    # Manually define segments
    segs = NULL
    segs = rbind(segs, c(arm_data$CN[1], arm_data$LOH[1], 1))
    for (r in 2:nrow(arm_data)){
      if (abs(arm_data$CN[r] - segs[nrow(segs),1]) > 0.05 | arm_data$LOH[r] != segs[nrow(segs),2]) {
        segs = rbind(segs, c(arm_data$CN[r], arm_data$LOH[r], r))
      }
    }
    segs = rbind(segs, c(arm_data$CN[r], arm_data$LOH[r], r))
    
    # Iterate through manual segments
    for (s in 1:(nrow(segs)-1)) {
      ss = segmentation(arm_data$BAF[segs[s,3]:segs[s+1,3]], n=n, g_size=g_size/2)
      start = which(baf$ARM == arm)[1]-1 + segs[s,3]-1 
      for (i in 1:nrow(ss)) {
        if ((ss[i,2]-ss[i,1] > s_limit) | (ss[i,2]-ss[i,1]) > length(which(baf$ARM == arm))*0.95) {
          segments((ss[i,1]+start)*scalefactor, ss[i,3:4]/100+ym, 
                   (ss[i,2]+start)*scalefactor, lwd=1.4, col = "#333333", xpd=T)
        }
      }
    }
    
  }
  
  dev.off()
  system(paste0("open ", result_file))
  
}
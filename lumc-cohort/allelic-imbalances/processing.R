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
lumc_data = readxl::read_xlsx("data/Supplementary Tables.xlsx", sheet=3)

# Read raw files
files = dir("lumc-cohort/allelic-imbalances/raw/")

# Iterate through samples
for (s_id in which(lumc_data$GS_ID %in% files[7])) {
  
  # Load sample/ID information
  sample = lumc_data$GS_ID[s_id]
  LUMC_id = lumc_data$ID[s_id]
  
  # Load BAF information
  baf = read.table(paste0("lumc-cohort/allelic-imbalances/raw/", sample, "/", sample, "_baf.tsv"))
  
  # Exclude X chromosome
  baf = baf[which(baf$CHR != "X"),]
  #baf = baf[which(baf$DP>quantile(baf$DP,0.4) & baf$DP<quantile(baf$DP,1)),]
  
  # Figure without segmentation
  result_file = paste0("lumc-cohort/allelic-imbalances/res/", LUMC_id, ".png")
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
    ymax = 1
    yat = c(1,2:(2+ymax))/6-(5+ymax)/6
    #segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
    axis(side = 2,at=chrmove+0.1,las=2,labels=c("Chromosome"),col.ticks = "#b1b1b1", 
         col = "#b1b1b1", col.axis="#333333", lwd=1.4,cex=1.1,cex.axis=1.1)
 
    # Define scalefactor (depending on total number of BAFs)
    scalefactor = 1/nrow(baf)*20
    
       # Dotted lines between each chromosomal arm
    for (arm in unique(baf$ARM)) {
      min = min(which(baf$ARM == arm))
      segments(min*scalefactor,0+ym,min*scalefactor,1+ym,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
      #segments(min*scalefactor,min(yat),min*scalefactor,-0.5,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
      #segments(min*scalefactor,min(yat)-0.15,min*scalefactor,min(yat)-0.15-1/6,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
    }
    
    # Solid lines between each chromosome + labels
    chromosome_data = table(baf$CHR)
    chromosome_data = chromosome_data[match(1:22, names(chromosome_data))]
    pos = 0.5
    for (chr in 1:length(chromosome_data)) {
      chromosome = names(chromosome_data)[chr]
      pos_new = pos + chromosome_data[chr]  
      
      segments(pos_new*scalefactor,0+ym,pos_new*scalefactor,1+ym,col = "#b1b1b1",lwd=1.4,xpd=T)
      #segments(pos_new*scalefactor,min(yat),pos_new*scalefactor,-.5,col = "#b1b1b1",lwd=1.4,xpd=T)
      #segments(pos_new*scalefactor,min(yat)-0.15,pos_new*scalefactor,min(yat)-0.15-1/6,col = "#b1b1b1",lwd=1.4,xpd=T)
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
    #segments(0,min(yat)-0.15,0,min(yat)-0.15-1/6,lwd=1.4,col="#b1b1b1",xpd=T)
    #segments(0,2/6-(5+ymax)/6,20,col="#b1b1b1",lwd=1.4)
    #segments(0,chrmove,20,col="#eeeeee",lwd=1.4)
    #segments(0,chrmove+0.2,20,col="#eeeeee",lwd=1.4)
    segments(0,chrmove+0.2,0,chrmove,col="#b1b1b1",lwd=1.4)
    # Title
    y=1.5
    text(-0.2, y, col="#333333", labels=LUMC_id, xpd=T, cex=1.1, font=2, pos=4) 
    
    # Labels
    x = 16.4
    text(x+0.3,y, cex=1.1,col="#333333", labels="BAF segmentation",xpd=T, pos=4)
    xx = c(0.08,0.18,0.00,0.03,0.23,0.13,0.09,0.23,0.12,0.12,0.26,0.09,0.18,0.10,0.14,0.20,0.14,0.10,0.03,0.26,0.21,0.26,0.16,0.13,0.13,0.17,0.21,0.07,0.12,0.19)
    yy = c(1.51,1.51,1.52,1.52,1.51,1.46,1.48,1.48,1.52,1.51,1.47,1.48,1.50,1.51,1.54,1.50,1.49,1.50,1.47,1.49,1.51,1.52,1.47,1.47,1.51,1.47,1.51,1.55,1.51,1.51)
    points(x+xx+0.01,yy,pch = 16,
           cex = 0.3,col=rgb(0,0,0,0.2),xpd=T)
    segments(x,y,x+0.3, lwd=1.4, col = "#333333", xpd=T)
    
  }
  
  dev.off()
  system(paste0("open ", result_file))
  
  # Figure with automated segmentation
  result_file = paste0("lumc-cohort/allelic-imbalances/res/", LUMC_id, "-segmented.png")
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
    ymax = 1
    yat = c(1,2:(2+ymax))/6-(5+ymax)/6
    #segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
    axis(side = 2,at=chrmove+0.1,las=2,labels=c("Chromosome"),col.ticks = "#b1b1b1", 
         col = "#b1b1b1", col.axis="#333333", lwd=1.4,cex=1.1,cex.axis=1.1)
    
    # Define scalefactor (depending on total number of BAFs)
    scalefactor = 1/nrow(baf)*20
    
    # Dotted lines between each chromosomal arm
    for (arm in unique(baf$ARM)) {
      min = min(which(baf$ARM == arm))
      segments(min*scalefactor,0+ym,min*scalefactor,1+ym,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
      #segments(min*scalefactor,min(yat),min*scalefactor,-0.5,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
      #segments(min*scalefactor,min(yat)-0.15,min*scalefactor,min(yat)-0.15-1/6,col = "#b1b1b1",lwd=1.4,xpd=T,lty=3)
    }
    
    # Solid lines between each chromosome + labels
    chromosome_data = table(baf$CHR)
    chromosome_data = chromosome_data[match(1:22, names(chromosome_data))]
    pos = 0.5
    for (chr in 1:length(chromosome_data)) {
      chromosome = names(chromosome_data)[chr]
      pos_new = pos + chromosome_data[chr]  
      
      segments(pos_new*scalefactor,0+ym,pos_new*scalefactor,1+ym,col = "#b1b1b1",lwd=1.4,xpd=T)
      #segments(pos_new*scalefactor,min(yat),pos_new*scalefactor,-.5,col = "#b1b1b1",lwd=1.4,xpd=T)
      #segments(pos_new*scalefactor,min(yat)-0.15,pos_new*scalefactor,min(yat)-0.15-1/6,col = "#b1b1b1",lwd=1.4,xpd=T)
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
    #segments(0,min(yat)-0.15,0,min(yat)-0.15-1/6,lwd=1.4,col="#b1b1b1",xpd=T)
    #segments(0,2/6-(5+ymax)/6,20,col="#b1b1b1",lwd=1.4)
    #segments(0,chrmove,20,col="#eeeeee",lwd=1.4)
    #segments(0,chrmove+0.2,20,col="#eeeeee",lwd=1.4)
    segments(0,chrmove+0.2,0,chrmove,col="#b1b1b1",lwd=1.4)
    # Title
    y=1.5
    text(-0.2, y, col="#333333", labels=LUMC_id, xpd=T, cex=1.1, font=2, pos=4) 
    
    # Labels
    x = 16.4
    text(x+0.3,y, cex=1.1,col="#333333", labels="BAF segmentation",xpd=T, pos=4)
    xx = c(0.08,0.18,0.00,0.03,0.23,0.13,0.09,0.23,0.12,0.12,0.26,0.09,0.18,0.10,0.14,0.20,0.14,0.10,0.03,0.26,0.21,0.26,0.16,0.13,0.13,0.17,0.21,0.07,0.12,0.19)
    yy = c(1.51,1.51,1.52,1.52,1.51,1.46,1.48,1.48,1.52,1.51,1.47,1.48,1.50,1.51,1.54,1.50,1.49,1.50,1.47,1.49,1.51,1.52,1.47,1.47,1.51,1.47,1.51,1.55,1.51,1.51)
    points(x+xx+0.01,yy,pch = 16,
           cex = 0.3,col=rgb(0,0,0,0.2),xpd=T)
    segments(x,y,x+0.3, lwd=1.4, col = "#333333", xpd=T)
    
  }
  for (arm in unique(baf$ARM)) {
    if (length(which(baf$ARM == arm))>5) {
    ss = segmentation(baf$BAF[which(baf$ARM == arm)], n=n*2, g_size=g_size*2)
    start = which(baf$ARM == arm)[1]-1
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

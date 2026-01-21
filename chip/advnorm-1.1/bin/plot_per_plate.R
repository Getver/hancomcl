# command line: R --no-save < code/plot_per_plate.R

############################
## Copyright (c) 2017 Thermo Fisher Scientific.  All Rights Reserved. This source
## code and the methods embodied within it represent the confidential
## intellectual property of Thermo Fisher Scientific, and this intellectual
## property is communicated on a confidential basis to the recipient.  Neither
## this code nor the methods may not be disclosed to any other party in any
## form without the express written permission of Thermo Fisher Scientific.  This
## copyright and notice must be preserved in any derivative work.
############################

options <- commandArgs(trailingOnly = TRUE)

if (length(options) != 6) {
  cat(paste("Usage: platemap.txt AxiomGT1.calls.txt AxiomGT1.summary.txt probesets_for_plate_normalization.ids output_dir [before|after]\n"))
  stop();
}

plate_map_file <- options[1]
callfile <- options[2]
sumfile <- options[3]
probeset_id_file <- options[4]
output_dir <- options[5]
status <- options[6]

############ read in files for plate information
##### d<-read.table("/nfs/falcon-usr/jinliu/Kaiser/data/CEL_header_info.txt",header=T,sep="\t")      # 480
d<-read.table(plate_map_file,header=T,sep="\t", colClasses = "character")      # 480
# summary(as.factor(d[,"plate_ID"]))
# summary(as.factor(d[,"plate"]))
cels <- d

#####################################
#####################################
# define other input fils
# pids: one col (col name: probeset_id), with the list of probesets ids that want to plot

##### pids <- read.table("/nfs/falcon-usr/jinliu/Kaiser/plot_slides/GWAS/pid.txt", as.is=T, header=T)
pids <- read.table(probeset_id_file, as.is=T, header=T)
pids<-pids[,1]

################## change these files: make a small file only contains the probesets that will be plot
##### indir <- "/nfs/falcon-usr/jinliu/Kaiser/genotyping_norm_v3/"
##### sumfile  <- paste(indir, "KP_GO_v2_058_norm.txt", sep="")
##### sumfile  <- '/nfs/falcon/usr/mingwu/advanced_axiom/plate_normalization/fas_training_201701/normalized_results/1_summary_intensities/1_normalized_summary.txt'
##### callfile <- paste(indir, "KP_GO_v2_058_SSP/AxiomGT1.calls.txt", sep="")
##### callfile <- '/nfs/falcon/usr/mingwu/advanced_axiom/plate_normalization/fas_training_201701/normalized_results/2_normalized_apt_results/normalized_results/AxiomGT1.calls.txt'
######################################################################################################

# cels <- d[,c("name","plate_ID")]
##### colnames(cels)<-c("cel_files","plate")
##### cels[,"cel_files"]<-paste(as.character(cels[,"cel_files"]),".CEL",sep="")
cels$plate <- as.factor(cels$plate)

#sam <- names(read.delim(pipe(paste('grep -v "^#" ', sumfile, ' | head -1', sep='')), header=T, check.names=F))[-1] 
sam <- read.table(callfile,header=F,nrow=1,sep="\t") 
sam<-as.character(as.matrix(sam[1,]))
sam<-sam[2:length(sam)]        # 1899 in the call file

rownames(cels)<-as.character(cels[,"cel_files"])
cels<-cels[as.character(sam),]

col <- rep(1, times=length(sam))



############ define functions
plot1 <- function(plate, p, pid, sumfile, inten) {
  a <- log2(as.numeric(inten[1,]))
  b <- log2(as.numeric(inten[2,]))
  M <- a - b
  A <- (a + b) * 0.5
  xlim <- max(abs(M))
  ymin <-min(A)
  ymax <-max(A)
  plot(M[plate==0], A[plate==0],
       col=8, cex=1.2, pch=16, 
       xlim=c(-xlim, xlim), ylim=c(ymin,ymax),xlab='', ylab='', axes=F, main=p)
  points(M[plate==1], A[plate==1],
       col=6, cex=1.2, pch=16, 
       xlim=c(-xlim, xlim), xlab='', ylab='', main='')
}

plotByCall <- function(pid, sam, sumfile, callfile, inten) {
  calls <- read.delim(pipe(paste('grep ', pid, ' ', callfile, sep='')),
  	   		header=F, colC=c('NULL', rep('numeric', length(sam))))

  colors <- rep(8, dim(calls)[[2]])
  colors[calls[1,] == 0] <- "red"
  colors[calls[1,] == 1] <- "gold"
  colors[calls[1,] == 2] <- "blue"

  a <- log2(as.numeric(inten[1,]))
  b <- log2(as.numeric(inten[2,]))
  M <- a - b
  A <- (a + b) * 0.5
  xlim <- max(abs(M))
  plot(M, A, col=colors, xlim=c(-xlim, xlim), xlab='', ylab='', axes=F,
       main=pid, cex=1.2, pch=16)
}

############## main function

num.plates <- length(unique(as.factor(d[,"plate"])))
num.rows <- ceiling(num.plates / 6)


for (j in 1:length(pids)) {
  pid <- pids[j]
  file <- file.path(output_dir, paste(pid, status, "png", sep="."))

  # png(file, width=1200, height=800)
  png(file, width=1200, height=(200 * num.rows))
  # par(mfrow=c(4,6), mar=c(1,1,1,1))
  par(mfrow=c(num.rows,6), mar=c(1,1,1,1))
  print(pid)

  inten <- read.delim(pipe(
     paste('grep ', pid, ' ', sumfile, sep='')),
     header=F, colC=c('NULL', rep('numeric', length(sam))))


  for (p in levels(as.factor(as.character(cels$plate)))) {
    plate <- rep(0, times=length(sam))
    for (i in 1:length(sam)) {
      if (cels$plate[cels$cel_files==sam[i]] == p) {
        plate[i] <- 1
      }
    }
    plot1(plate, p, pid, sumfile,inten)
  } 

  #plotByCall(pid, sam, sumfile, callfile, inten)


  dev.off()
}

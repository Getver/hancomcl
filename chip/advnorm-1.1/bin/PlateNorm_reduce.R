#!/usr/bin/env Rscript

############################
## Copyright (c) 2017 Thermo Fisher Scientific.  All Rights Reserved. This source
## code and the methods embodied within it represent the confidential
## intellectual property of Thermo Fisher Scientific, and this intellectual
## property is communicated on a confidential basis to the recipient.  Neither
## this code nor the methods may not be disclosed to any other party in any
## form without the express written permission of Thermo Fisher Scientific.  This
## copyright and notice must be preserved in any derivative work.
############################

###############################  
# processing local files
###############################

# pidlist <-"platenorm_pidlist.txt"
# plate_info <- "platemap.txt"
# sam_order <- "sam_order.txt"      # one file for each batch
# sample_mask <- "sample_mask.txt"

args <- commandArgs(trailingOnly = TRUE)
if ((length(args) > 4) || (length(args) < 3)) {
  stop("\n\nArgument 2: probeset_ids_for_normalization_file\nArgument 2: platemap_file (required)\nArgument 3: sample_order_file (required)\nArgument 4: sample_mask_file (optional)")
}

pidlist <- args[1]
plate_info <- args[2]
sam_order <- args[3]

mask.indexes <- c()

if (length(args) == 4) {
  sample_mask <- args[4]
  mask <- read.table(sample_mask, header=T)
  mask.indexes <- as.integer(mask[,1]) + 1
}

sam_min <- 8

sam <- read.table(sam_order,header=F)
sam <- as.matrix(sam)
sam <-as.character(sam[,1])

mask.flag=FALSE
if (length(mask.indexes > 0)) {
  mask.indexes = mask.indexes * -1
  sam <- sam[mask.indexes]
  mask.flag=TRUE
}

# read plate info from a file
cels<-read.table(plate_info,header=T,sep="\t", colClasses = c("character", "character"))      # 480
cels<-as.matrix(cels)
colnames(cels)<-c("cel_files","plate")
# make sure cels & summary files have the same sample order
rownames(cels)<-cels[,"cel_files"]
cels <- cels[as.character(sam),]
plate <- as.factor(cels[,"plate"])

# set small plates to NA (therefore no contribution in the analysis)
psize<-summary(as.factor(cels[,"plate"]))
pkeep<-names(psize)[psize>sam_min]
plate<-ifelse(as.character(plate) %in% pkeep, as.character(plate), NA)
plate <- as.factor(plate)

# read in the pid list and only do normalization for the probesets in this list
ref<-read.table(pidlist,header=F)
ref_A<-paste(as.character(ref[,1]),"-A",sep="")
ref_B<-paste(as.character(ref[,1]),"-B",sep="")
ref_AB<-as.character(c(as.character(ref_A),as.character(ref_B)))


###############################
# define the function
###############################
norm <- function(value,plate){
  inten <- value
  pid <-inten[1]

  inten <-inten[2:length(inten)]
  if (isTRUE(mask.flag)) {
    inten <- inten[mask.indexes]
  }

  inten <- as.integer(inten)
            
  r <- lm(log(inten)~ plate) 
  inten_norm <- exp(r$coefficients[1] + r$residuals)
  names(inten)<-c(1:length(inten))
  inten[names(inten_norm)]<-inten_norm
                   
  inten_norm <- round(inten,5)
  inten_norm <-c(as.character(pid),as.character(inten_norm))
  result <-paste(inten_norm,collapse="\t")
  
  # return a single value or a vector 
  return(result)
}


###############################
# stdin the output from map
###############################

# reducer
current.key <- NA
current.val <- 0.0

conn <- file("stdin",open="r")               
while(length(next.line <- readLines(conn, n=1)) > 0){
  split.line <- strsplit(next.line, "\t")
  key <- split.line[[1]][1]
  # if(key %in% as.character(ref_AB)){
  if(key %in% ref_AB){
    val <- split.line[[1]]
  
    if(is.na(current.key)){
      current.key <- key
      current.val <- val
    }
    else{
      if(current.key == key){ 
      }
      else{   
        inten_norm <-norm(current.val,plate)           
        write(inten_norm,stdout())
             
        current.key <- key
        current.val <- val
      }
    }
  }
}

if (!(is.na(current.key))) {
  inten_norm <-norm(current.val,plate)           
  write(inten_norm,stdout())
}

close(conn)

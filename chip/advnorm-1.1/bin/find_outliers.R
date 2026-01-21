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
## processing local files
################################

platedet.sam_min <- 8

# command-line argument handling

args <- commandArgs(trailingOnly = TRUE)
if ((length(args) > 3) || (length(args) < 2)) {
  stop("\n\nArgument 1: platemap_file (required)\nArgument 2: sample_order_file (required)\nArgument 3: sample_mask_file (optional)")
}

platedet.plate_info <- args[1]
platedet.sam_order <- args[2]

platedet.mask.indexes <- c()

if (length(args) == 3) {
  platedet.sample_mask <- args[3]
  platedet.mask <- read.table(platedet.sample_mask, header=T)
  platedet.mask.indexes <- as.integer(platedet.mask[,1]) + 1
}

platedet.sam <- read.table(platedet.sam_order,header=F)
platedet.sam <- as.matrix(platedet.sam)
platedet.sam <-as.character(platedet.sam[,1])

if (length(platedet.mask.indexes > 0)) {
  platedet.mask.indexes = platedet.mask.indexes * -1
  platedet.sam <- platedet.sam[platedet.mask.indexes]
}

platedet.sam.length <- length(platedet.sam)

# read plate info from a file
platedet.cels<-read.table(platedet.plate_info,header=T,sep="\t", colClasses = c("character", "character"))      # 480
platedet.cels<-as.matrix(platedet.cels)
colnames(platedet.cels)<-c("cel_files","plate")
# make sure cels & summary files have the same sample order
rownames(platedet.cels)<-platedet.cels[,"cel_files"]
platedet.cels <- platedet.cels[as.character(platedet.sam),]
platedet.plate <- as.factor(platedet.cels[,"plate"])

platedet.psize<-summary(as.factor(platedet.cels[,"plate"]))
platedet.pkeep<-names(platedet.psize)[platedet.psize>platedet.sam_min]
platedet.plate<-ifelse(as.character(platedet.plate) %in% platedet.pkeep, as.character(platedet.plate), NA)
platedet.plate <- as.factor(platedet.plate)

header <- paste(c('probeset_id', 'outlier_status', 'outlier_plate_ids'), collapse="\t")
write(header, stdout())

conn <- file("stdin",open="r")
while(length(calls.s <- readLines(conn, n=1)) > 0){

  calls.s <- strsplit(calls.s, "\t", fixed=TRUE)[[1]]

  # pid <- top.probeset_id
  pid <- calls.s[1]
  call <- calls.s[2:length(calls.s)]

  # we need to exclude the calls for samples that we wish to exclude
  # since this version isn't working the same way as the version
  # we need for Hadoop processing

  if (length(platedet.mask.indexes > 0)) {
    call <- call[platedet.mask.indexes]
  }

  if (length(call) != platedet.sam.length) {
    stop("Calls input does not have expected number of samples")
  }

  call <- as.integer(call)

  tag<- !(call==-1)    

  # reads from a global!
  t<-table(platedet.plate[tag],call[tag])
  t1<-prop.table(t,margin=2)*100             
  t2<-prop.table(t,margin=1)*100   
  ord<-rownames(t)

  # initialize
  # result="NO PLATE DETECTION RESULTS"
  result=NA

  if(ncol(t)==1){
    result <- rep(0,nrow(t))
    # result <-c(as.character(pid),as.character(result))
    # result <-paste(result,collapse="\t")
  } else if(ncol(t)==2){
    tag_col<-colSums(t)    
    tag_col <- (tag_col==min(tag_col))  
    tag_col <- ifelse(duplicated(tag_col),FALSE,tag_col)

    if(sum(t[,tag_col])>=10){  # have to have >=5 points in the smaller cluster
      t2<-t2[,tag_col]
      r_c<-ifelse(t2>70,1,0)
      t1<-t1[,tag_col]
      r_p<-ifelse(t1>10,1,0)
      r_add<-ifelse(t1>60,1,0)
      r<- (r_c & r_p)
      r<- (r | r_add)
      result<-ifelse(r=="TRUE",1,ifelse(r=="FALSE",0,as.character(r)))
    } else {
      result <- rep(0,nrow(t))      
    } 
        
    # result <-c(as.character(pid),as.character(result))
    # result <-paste(result,collapse="\t")

  } else if(ncol(t)==3){
    tag_col<-colSums(t)    
    minc <- (tag_col==min(tag_col))
    minc <- ifelse(duplicated(minc),FALSE,minc)
    maxc <- (tag_col==max(tag_col))
    maxc <- ifelse(minc==TRUE,FALSE,maxc)
    maxc <- ifelse(duplicated(maxc),FALSE,maxc)
    midc <- !(minc | maxc)
         
    ta <- t[, (minc | maxc)]
    t1<-prop.table(ta,margin=2)*100      
    t1<-t1[as.character(ord),]    
    t2<-prop.table(ta,margin=1)*100      
    t2<-t2[as.character(ord),]
    t1<-ifelse(is.na(t1),0,t1)
    t2<-ifelse(is.na(t2),0,t2)

    tag_col<-colSums(ta) 
    tag_col <- (tag_col==min(tag_col))  
    tag_col <- ifelse(duplicated(tag_col),FALSE,tag_col)
    if(sum(ta[,tag_col])>=10){  # have to have >=5 points in the smaller cluster
      t2<-t2[,tag_col]
      r_c<-ifelse(t2>70,1,0)
      t1<-t1[,tag_col]
      r_p<-ifelse(t1>10,1,0)
      r_add<-ifelse(t1>60,1,0)
      r<- (r_c & r_p)
      r_a<- (r | r_add)
    } else {
      r_a <- rep(FALSE,nrow(t))      
    } 

    tb <- t[, (midc | maxc)]
    t1<-prop.table(tb,margin=2)*100      
    t1<-t1[as.character(ord),]        
    t2<-prop.table(tb,margin=1)*100      
    t2<-t2[as.character(ord),]
    t1<-ifelse(is.na(t1),0,t1)
    t2<-ifelse(is.na(t2),0,t2)

    tag_col<-colSums(tb) 
    tag_col <- (tag_col==min(tag_col))  
    tag_col <- ifelse(duplicated(tag_col),FALSE,tag_col)
    if(sum(ta[,tag_col])>=10){  # have to have >=5 points in the smaller cluster
      t2<-t2[,tag_col]
      r_c<-ifelse(t2>70,1,0)
      t1<-t1[,tag_col]
      r_p<-ifelse(t1>10,1,0)
      r_add<-ifelse(t1>60,1,0)
      r<- (r_c & r_p)
      r_b<- (r | r_add)
    } else {
      r_b <- rep(FALSE,nrow(t))      
    } 

    r <- r_a | r_b

    result<-ifelse(r=="TRUE",1,ifelse(r=="FALSE",0,as.character(r)))
    # result <-c(as.character(pid),as.character(result))
    # result <-paste(result,collapse="\t")

  }

  # write a single value or a vector 
  # print(class(result))
  # if (is.na(result)) {
  if (anyNA(result)) {
    result <- paste(c(as.character(pid), "NO PLATE DETECTION RESULTS"), collapse="\t")
  } else {
    status <- "NO OUTLIER PLATE DETECTED"
    outlier.plate.ids <- ''

    if (sum(as.integer(result)) > 0) {
      status <- "OUTLIER PLATE DETECTED"
      # write(status, stdout())
      # print(names(r))
      # print(r)
      outlier.plate.ids <- paste(c(names(r)[r]), collapse=",")
      # print(outlier.plate.ids)

    }

    result <- paste(c(as.character(pid), status, outlier.plate.ids, as.character(result)), collapse="\t")

  }

  write(result, stdout())

}

close(conn)

library(RCurl)
library(readr)

URL <- url("https://raw.githubusercontent.com/AbrahamProject/MZ_ERROR/main/RETENTION_TIME.csv")
standard <- read_csv(URL,show_col_types = FALSE)

standard <- data.frame(standard)




mass_error <- function(data,mz,mz_delta,rt_delta,SAVE=F,adduct){
  if(SAVE==T){  

  p <- mz
  standard <-standard
  c <- adduct
  mz <- unlist(c(p,c))
  mzr <- c(standard[mz[1],mz[2]]-mz_delta,standard[mz[1],mz[2]]+mz_delta)
  rtr <- rt_delta
  
  ms1 <- which(msLevel(data) == 1)
  rtsel <- rtime(data)[ms1] > rtr[1] & rtime(data)[ms1]  < rtr[2]
  
  
  M3 <- MSmap(data, ms1[rtsel], mzr[1], mzr[2], .005)
  
  
  datatable <- M3@map
  rownames(datatable) <- M3@rt
  colnames(datatable) <- M3@mz
  
  
  maximum<-matrix(ncol=2,nrow=ncol(datatable))
  maximum[,1] <- colnames(datatable)
  colnames(maximum) <- c("mz","#")
  for(i in 1:ncol(datatable)){
    
    if(all(datatable[,i]==0) ){
      maximum[i,2] <- NA
    }else{
      index <- which(datatable[,i]>0)
      maximum[i,2] <- length(index) 
    }
  }
  
  
  mz_error <- matrix(ncol=2,nrow=nrow(maximum))
  mz_error[,1] <- as.numeric(maximum[,1])-as.numeric(standard[mz[1],mz[2]]) 
  mz_error[,2] <- maximum[,2]
  
  
  mz_error <- mz_error[!is.na(mz_error[,2]),]
  
  return(mz_error)
  
  
  
}else{
  
  p <- mz
  standard <-standard
  c <- adduct
  mz <- unlist(c(p,c))
  mzr <- c(standard[mz[1],mz[2]]-mz_delta,standard[mz[1],mz[2]]+mz_delta)
  rtr <- rt_delta
  
  ms1 <- which(msLevel(data) == 1)
  rtsel <- rtime(data)[ms1] > rtr[1] & rtime(data)[ms1]  < rtr[2]
  
  
  M3 <- MSmap(data, ms1[rtsel], mzr[1], mzr[2], .005)
  
  
  
  datatable <- M3@map
  rownames(datatable) <- M3@rt
  colnames(datatable) <- M3@mz
  
  
  
  chrom_int <- data.frame(rt=rownames(datatable),int=rowSums(datatable))
  
  pal <- colorRamp(c("blue", "green", "orange", "red"))
  
  x<-chrom_int$int
  if(sum(x)!=0){
    if( length(which(x>0))>1 ){
  col1 <- rgb(pal((x - min(x)) / diff(range(x))), max=255)
  
  chrom_mz <-data.frame(rt=rep(rownames(datatable),each=ncol(datatable)),mz=rep(colnames(datatable),nrow(datatable)),int=c(t(datatable)))
  
  s<-chrom_mz[which(chrom_mz$int!=0),] 
  
  x<-s$int
  col2 <- rgb(pal((x - min(x)) / diff(range(x))), max=255)
  
  
  
  layout(matrix(c(1,2,3,3),nrow=2),width=c(5,5))
  par(mar=c(0,6,3,3))
  plot(chrom_int,col=col1,pch=16,cex=1.2,xlab="",xaxt="n",xlim=c(rtr[1],rtr[2]),main=paste0("mz accuracy: ",standard[mz[1],mz[2]]))
  par(mar=c(4.5,6,0,3))
  plot(s$rt,s$mz,col=col2,pch=16,las=1,xlim=c(rtr[1],rtr[2]),xlab="retention time [s]",ylab="mz",cex=1.2)
  
  
  
  
  ##########################################################################################################################
  
  maximum <- 0
  
  for(i in 6:9){
    mz <- p
    mzr <- c(standard[p,(i)]-mz_delta,standard[p,(i)]+mz_delta)
    ms1 <- which(msLevel(data) == 1)
    rtsel <- rtime(data)[ms1] > rtr[1] & rtime(data)[ms1]  < rtr[2]
    
    
    M3 <- MSmap(data, ms1[rtsel], mzr[1], mzr[2], .005)
    
    
    
    datatable <- M3@map
    rownames(datatable) <- M3@rt
    colnames(datatable) <- M3@mz
    
    
    
    chrom_int <- data.frame(rt=rownames(datatable),int=rowSums(datatable))
    maximum2 <-max(chrom_int$int)
    if(maximum2>maximum){
      maximum <- maximum2
    }
  }
  
  
  
  
  
  mzr <- c(standard[p,6]-mz_delta,standard[p,6]+mz_delta)
  ms1 <- which(msLevel(data) == 1)
  rtsel <- rtime(data)[ms1] > rtr[1] & rtime(data)[ms1]  < rtr[2]
  
  
  M3 <- MSmap(data, ms1[rtsel], mzr[1], mzr[2], .005)
  
  
  
  datatable <- M3@map
  rownames(datatable) <- M3@rt
  colnames(datatable) <- M3@mz
  
  
  
  chrom_int <- data.frame(rt=rownames(datatable),int=rowSums(datatable))
  
  color <- c("#000000","#D0CECE","#5B9BD5","#70AD47")
  par(mar=c(4,4,3,3))
  
  
  plot(chrom_int,xlab = "RT",ylab = "intenisty",type="l", col = color[1], ylim=c(0,1.05*maximum),lwd=2)
  
  
  
  for(i in 7:9){
    mzr <- c(standard[p,i]-mz_delta,standard[p,i]+mz_delta)
    ms1 <- which(msLevel(data) == 1)
    rtsel <- rtime(data)[ms1] > rtr[1] & rtime(data)[ms1]  < rtr[2]
    
    
    M3 <- MSmap(data, ms1[rtsel], mzr[1], mzr[2], .005)
    
    
    
    datatable <- M3@map
    rownames(datatable) <- M3@rt
    colnames(datatable) <- M3@mz
    
    
    
    chrom_int <- data.frame(rt=rownames(datatable),int=rowSums(datatable))
    
    
    points(chrom_int[,1],chrom_int[,2],col=color[i-5],type = "l",lwd=2)
  }
  
  legend("topright",legend = colnames(standard)[6:9], col=color,lwd=2,bty="n",cex=1.2 )
  }}
}
}
  




  


###Axoquant2.0

###Please replace the example text in quotes with the path to your experiment folder.
experiment.folder<-("/users/Aaron/Desktop/test")
source("http://bioconductor.org/biocLite.R")

###You will be asked to download Bioconductor by executing the next line. You will require an internet connection. 
###Please allow it to download before proceeding.

biocLite("EBImage")                     


library("EBImage")
setwd(experiment.folder)

Path <- experiment.folder
output_directory <- experiment.folder
colour <- "blue"

Image.process <-function(Path, output_directory, colour,threshold=NA){   ###To use a custom threshold value, replace NA with your value. Otherwise, each image will be auto-thresholded
 
  #Set up dataframe to store final data
  Finaloutput<-data.frame(1:76)                                
  Finalindiv<-data.frame(1:80)
  
  #Get folders for each embryo from   
  Directory<-list.dirs(Path,recursive=TRUE)
  Directory<-as.data.frame(Directory)
  Directory<-Directory[2:nrow(Directory),]
  Directory<-as.data.frame(Directory)

  for(a in 1:nrow(Directory)){
    
    setwd(paste(Directory[a,1], sep=""))                    
  
    files<-list.files(paste(Directory[a,1], sep=""))        
  
    #Init dataframe to store data for plot
    Plotdata<-data.frame(1:80)                             
    #Init dataframe to store Intensity/Distance data
    Data<-data.frame(1:10000)
    #Load an image and cycle through
    for(b in 1:length(files)) {
      iteration <- readImage(files[b])
      print(b,quote=F)
    
      size<-dim(iteration)
      widthscalar<-size[1]/100
      heightscalar<-size[2]/100
      #Resize
      resized<-resize(iteration,100,100)
      #Extract green channel
      matrix<-resized[1:100,1:100,2]
      
      if(is.na(threshold)){                              
        thresholded<-matrix> round(mean(iteration)+1*(sd(iteration)),4)
      }else{
        thresholded<-matrix>threshold
      }
    
      colorMode(thresholded)<-Grayscale
      #Get size, first number is widthm second is length
    
      #Calculated intensity of each corner
      upleft<-mean(thresholded[1:10,1:10])
      upright<-mean(thresholded[90:100,1:10])
      lowleft<-mean(thresholded[1:10, 90:100])
      lowright<-mean(thresholded[90:100,90:100])
      corners<-list(upleft,upright,lowleft,lowright)
      corners<-as.data.frame(corners)
      
      #Rotate so that bright corner is in upper left
      rotated<-thresholded
      if(max(corners)==upleft){}                                        
      if(max(corners)==upright){rotated<-rotate(rotated, 270)}
      if(max(corners)==lowleft){rotated<-rotate(rotated, 90)}
    
      if(max(corners)==lowright){rotated<-rotate(rotated, 180)}
    
      #Generate distance matrix, distance is in um
      dist.matrix<-matrix(data<-NA, nrow<-100, ncol=100)
    
      for(j in 1:100){
        for(i in 1:100){
          dist.matrix[i,j] = sqrt((i*widthscalar)^2+(j*heightscalar)^2)
        }
      }
    
      #Start calculating distance and intensity
      rotated<-as.data.frame(rotated)
      #init data.frame
      mydata<-data.frame(2,1:10000)
      #Start count
      k<-1
      #Loop through and fill out mydata matrix
      for(i in 1:100) {
        for(j in 1:100){
          mydata[k,1]<-rotated[i,j]
          mydata[k,2]<-dist.matrix[i,j]
          k<-k+1
        }
      }
    
      #Name data
      colnames(mydata) <- c("Intensity","Distance") 
      #Sort by distance
      ordered<-mydata[order(mydata$Distance),]
      #Save data in Data
      Data<-cbind(Data,ordered)
      #Bin
      binned<-tapply(ordered$Intensity, cut(ordered$Distance, seq(0, 4000, by=50)), mean)
      Plotdata<-cbind(Plotdata,binned)
    }
    
    Finaldata<-Plotdata[,2:(length(files)+1)]
    Indivtemp<-Finaldata
    colnames(Indivtemp)<-files
    Finalindiv<-cbind(Finalindiv,Indivtemp)  
  
    #Graph
    setwd(output_directory)  
    pdf(paste("Output","_",a,".pdf", sep=""))                        
                                                                     
    #Graph data from specified folder path
    Finaldata[is.na(Finaldata)] <- 0
    Finaldata<-Finaldata[5:80,]
    Temp<-Finaldata[,1]
    plot(Temp, type="l", main=paste("",Directory[a,1],".pdf", sep=""), xlab="Distance from soma in 50 pixel units", ylab="Axon density",ylim=c(0,1), xlim=c(0,80),col=colour)
    #Loop through to add lines
    for(i in 1:(dim(Finaldata)[2])){
      Temp<-Finaldata[,i]
      lines(Temp, col=colour)  
    }
    #Graph average
    Mean<-rowMeans(Finaldata, na.rm = FALSE, dims = 1)
    Mean<-as.data.frame(Mean)
    lines(Mean, col="black", lwd=5)
    legend('topright', c('DRG','Mean'), lty=1, lwd=c(1,5), col=c(colour,'black'), bty='o', cex=.75)
    dev.off()
  
    Finaldata$avg <- apply(Finaldata,1,mean,na.rm=TRUE) 
    Temp<-Finaldata[,ncol(Finaldata)]
    Temp<-as.data.frame(Temp)
    colnames(Temp)<-Directory[a,1]
  
    Finaloutput<-cbind(Finaloutput,Temp)
  }

  # remove 1st column in Finaloutput
  Finaloutput<-Finaloutput[,2:ncol(Finaloutput)]
  
  #Graph with only mean values
  pdf(file="Output_means.pdf")   
  plot(0,0, type="l", main="Mean values", xlab="Distance from soma in 50 micron units", ylab="Axon density",ylim=c(0,1), xlim=c(0,80),col="white")
  #Loop through to add lines
  colorsMean <- colorRampPalette(c("red","yellow","springgreen","royalblue"))(ncol(Finaloutput))   
  namesMean <- c()
  for(i in 1:ncol(Finaloutput)){
    Tempo<-Finaloutput[,i]
    lines(Tempo, col=colorsMean[i],lwd=4)
    namesMean[i] <- paste("folder_mean",i,sep="")
  }
  legend('topright', namesMean, lty=1, lwd=4, col=colorsMean, bty='o', cex=.75)
  dev.off()
  
  Finalindiv<-Finalindiv[5:80,2:ncol(Finalindiv)]
  Finalindiv[is.na(Finalindiv)] <- 0
  
  ## creates a list of 2 dataframes, one with the individual measurements and a second with the mean values
  FinalouputList <- list(Finalindiv,Finaloutput)
  names(FinalouputList) <- c("Individual_measurements","Mean_values")
  
  return(FinalouputList)
}

##Begin the image analysis
Parent <- Image.process(experiment.folder, experiment.folder,"blue",threshold=NA)

write.table(Parent$Individual_measurements,"Individual_measurements.csv", sep=",", row.names=FALSE)
write.table(Parent$Mean_values,"Embryo_Means.csv", sep=",", row.names=FALSE)

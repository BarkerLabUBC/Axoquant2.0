###Axoquant2.0 as of May 10 2018

###Please replace the example path in quotes with the path to your own experiment folder.
experiment.folder<-("/users/A/Desktop/Example")

source("http://bioconductor.org/biocLite.R")
###You will be asked to download Bioconductor by executing the next line. You will require an internet connection. 
###Please allow it to download before proceeding.

biocLite("EBImage")                     


library("EBImage")

setwd(experiment.folder)

Path <- experiment.folder
output_directory <- experiment.folder
colour <- "blue"

###### The 10000 one
Image.process <-function(Path, output_directory, colour,pixdim,threshold=NA){   ###To use a custom threshold value, replace NA with your value. Otherwise, each image will be auto-thresholded
  
  #Set up dataframe to store final data
  Finaloutput<-data.frame(1:196)                                
  Finalindiv<-data.frame(1:200)
  
  #Get folders for each embryo from   
  Directory<-list.dirs(Path,recursive=TRUE)
  Directory<-as.data.frame(Directory)
  Directory<-Directory[2:nrow(Directory),]
  Directory<-as.data.frame(Directory)
  
  for(a in 1:nrow(Directory)){
    ### Change back to [a,1]
    setwd(paste(Directory[a,1], sep=""))                    
    ### Change back to [a,1]
    files<-list.files(paste(Directory[a,1], sep=""))        
    
    #Init dataframe to store data for plot
    Plotdata<-data.frame(1:200)                             
    #Init dataframe to store Intensity/Distance data
    Data<-data.frame(1:100000)
    #Load an image and cycle through
    for(b in 1:length(files)) {
      ### Change back to [b]
      iteration <- readImage(files[b])
      print(b,quote=F)
      
      size<-dim(iteration)
      widthscalar<-size[2]/100
      heightscalar<-size[1]/100
      #Resize
      resized<-resize(iteration,100,100)
      #Extract green channel
      matrix<-resized[1:100,1:100,2]
      #### REmove next 2 lines
  
      
      if(is.na(threshold)){                              
        thresholded<-matrix> round(mean(iteration)+1.5*(sd(iteration)),4)
      }else{
        thresholded<-matrix>threshold
      }
      
      colorMode(thresholded)<-Grayscale
      #Get size, first number is width second is length
      
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
          dist.matrix[i,j] = sqrt((i*widthscalar*pixdim)^2+(j*heightscalar*pixdim)^2)
        }
      }
      
      #Start calculating distance and intensity
      rotated<-as.data.frame(rotated)
      #init data.frame
      mydata<-data.frame(0,1:10000)
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
      binned<-tapply(ordered$Intensity, cut(ordered$Distance, seq(0, 4000, by=20)), mean) 
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
    Finaldata<-Finaldata[5:200,]
    Temp<-Finaldata[,1]
    plot(Temp, type="l", main=paste("",Directory[a,1],".pdf", sep=""), xlab="Distance from soma (bin no.)", ylab="Axon density",ylim=c(0,1), xlim=c(0,200),col=colour)
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
  #pdf(file="Output_means.pdf")   
  #plot(0,0, type="l", main="Mean values", xlab="Distance from soma in 500 micrometer units", ylab="Axon density",ylim=c(0,1), xlim=c(0,80),col="white")
  #Loop through to add lines
  #colorsMean <- colorRampPalette(c("red","yellow","springgreen","royalblue"))(ncol(Finaloutput))   
  #namesMean <- c()
  #for(i in 1:ncol(Finaloutput)){
  #  Tempo<-Finaloutput[,i]
  #  lines(Tempo, col=colorsMean[i],lwd=4)
  #  namesMean[i] <- paste("folder_mean",i,sep="")
  #}
  #legend('topright', namesMean, lty=1, lwd=4, col=colorsMean, bty='o', cex=.75)
  #dev.off()
  
  #Finalindiv<-Finalindiv[5:200,2:ncol(Finalindiv)]
  #Finalindiv[is.na(Finalindiv)] <- 0
  
  ## creates a list of 2 dataframes, one with the individual measurements and a second with the mean values
  FinalouputList <- list(Finalindiv,Finaloutput)
  names(FinalouputList) <- c("Individual_measurements","Mean_values")
  
  return(FinalouputList)
}

##Begin the image analysis

Parent <- Image.process(experiment.folder, experiment.folder,"blue",3.87, threshold=NA) ##replace 3.87 with your own pixel to micrometer scale. 


write.table(Parent$Individual_measurements,"Individual_measurements.csv", sep=",", row.names=FALSE)
write.table(Parent$Mean_values,"Embryo_Means.csv", sep=",", row.names=FALSE)

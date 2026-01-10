# install.packages("ape")
# install.packages("stats")
# install.packages("utils")
# install.packages("ggplot2")
# install.packages("parallel")
# install.packages("pbapply")

library(ape)
library(stats)
library(utils)
library(ggplot2)
library(parallel)
library(pbapply)


fasta_to_table <- function(Seq_file){
  
  
  #Initializing the flag
  fasta_flag=0
  Header<-NULL
  Sequence<-NULL
  
  for (j in 1:nrow(Seq_file)){
    
    #This is setting up a flag so the first time we initialize the matrix and then the second time we rbind to the matrix for the fasta file format
    if (fasta_flag==0){
      if(grepl(">",Seq_file[j,1])==TRUE){
        Header<-as.vector(Seq_file[j,1])
        seq_concate<-""
        fasta_flag=1
      }
    }else{
      if(grepl(">",Seq_file[j,1])==TRUE){
        Sequence<-rbind(Sequence,seq_concate)
        Header<-rbind(Header,as.vector(Seq_file[j,1]))
        seq_concate<-""
      }else{
        seq_concate<-paste(seq_concate,Seq_file[j,1],sep="")
      }
    }
  }
  Sequence<-rbind(Sequence,seq_concate)
  
  two_column_fasta<-cbind(Header,Sequence)
  
  return(two_column_fasta)
  
}

AA_code = "invert"
AGCT_only = TRUE
data_folder = NULL
gen_outliers = TRUE
sp_outliers = TRUE
dist_model = "raw"
replicates = 0
replacement = TRUE
conf_level = 1
numCores = 1



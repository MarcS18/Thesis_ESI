##This R script calculates Q2 (predictivity) of successive principal components for a dataset
##INPUTS:
##INPUTS:
##Descriptors.csv - a .csv file containing calculated descriptors for a dataset, the descriptors must be named as in the global variable "Descs"
##OUTPUTS:
##Descriptors_Q2.csv - a .csv containing the Q2 values for successive principal components, first column is the principal component number and the second column is the Q2 value

##section 1: import libraries
library("pcaMethods")
library("funr")

wd <- getwd()
setwd(wd)
##section 2: define inputs and outputs
dir <- funr::get_script_path()##get path to directory .r script is in
Descriptors <- "scripts/ml/TSSig_PYR_full_trimmed.csv" ##name of file containing descriptors
Descriptors_Q2 <- "scripts/ml/Descriptors_TSSig_PYR_Q2.csv"##name of output file

##section 3: define method for getting Q2
##column names of descriptors
sdata <- data.frame(x = names(read.csv(file.path(wd,Descriptors))))

# sdata <- read.table(file.path(wd,Descriptors),             # Read only header of example data
#            head = TRUE,
#            nrows = 1,
#            sep = ",")[- 1, - 2]
descs <- sdata$x
descs <- head(descs, -2)
descs <- descs[-1]
#define method
get_Q2 <- function(dir,Descriptors,Descriptors_Q2){
  ##read in data
  data <- read.csv(file.path(dir,Descriptors))
  ##get just descriptors
  x <- data[descs]
  #run PCA analysis
  pcIr <- pca(x, nPcs=length(descs), cv = "q2", scale="vector")
  ##get Q2 values, Components below a certain threshold are excluded
  q2 <- Q2(pcIr, x)
  ##write output to file
  write.csv(q2,file.path(dir,Descriptors_Q2), row.names = FALSE)
}

##section 4: run method to get Q2
get_Q2(wd,Descriptors,Descriptors_Q2)
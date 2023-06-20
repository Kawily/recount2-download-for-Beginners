
###--------------------------------DON`T PANIC -----------------------------###


#This is the recount2 script

## Install recount from Bioconductor
install.packages("BiocManager")
BiocManager::install('recount')

## Browse the vignetets for a quick description of how to use the package
library('recount')
browseVignettes('recount')

## Download the RangedSummarizedExperiment object at the gene level for 
## study SRP045638
url <- download_study('SRP045638')

## View the url for the file by printing the object url
url

## Load the data
load(file.path('SRP045638', 'rse_gene.Rdata'))

## Scale counts
rse <- scale_counts(rse_gene)

#The next step is done in Excel not in R because Excel makes out of a complex
#csv a simple xl data, which kinda solved my format problem...
#So first create a csv
counts_brain <- rse@assays@data@listData[["counts"]]

write.csv(counts_brain,"C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit/counts_brain.csv" )
library(readr)
setwd("C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit")
rcounts <- read.csv("counts_brain.csv")

#Second open it on excel:Daten/Datenabrufen/AusDatei/ausCSV (Dann auf Laden NICHT TRANSFORMIEREN!!)
#Save it as xlsx data (e.g. Book1.xlsx)
#If it all worked we have a Dataformat that is easy to work with

library(readxl)
setwd("C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit")
rawcounts <- read_xlsx("Book1.xlsx")

#we have to change the Ensemble IDs from first row to row name
#And we have to get rid of the ENS...0000.12 version at the end
#And we have to change the format from tibble to dataframe

#Change tibble to data frame
str(rawcounts)
rawcounts <- as.data.frame(rawcounts)
str(rawcounts)

raw <- rawcounts[,-1]
rawnames <- rawcounts[1]

#we get rid of the ENS Version via splitting the column

library(dplyr)
library(tidyr)

rnamesversionfree <-pivot_longer(rawnames, cols = everything()) %>%
  separate(value, into = c('gene_id', 'beta'), sep = "\\.") %>%
  select(-name)

#Now we only have to make the rnamesversionfree to our row names
#some format changes are needed

rnames <- rnamesversionfree[,1]
rnamesdata <- as.data.frame(rnames)

rawnamesdataunique <-make.names(rnamesdata[,1], unique = TRUE)

row.names(raw) <- rawnamesdataunique

#be happy -> the raw data is now perfectly fine to work with RNAAGECalc functions


#############-----For Huge Data e.g. GTEX --------###########

url <- download_study('SRP012682')

## View the url for the file by printing the object url
url

## Load the data
load(file.path('SRP012682', 'rse_gene.Rdata'))

## Scale counts
rse <- scale_counts(rse_gene)


#So first create a csv...
counts_brain <- rse@assays@data@listData[["counts"]]
setwd("C:/Users/Geologie5.ALPECON/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit")
write.csv(counts_brain,"(GTEX_data1.csv" )

#This is when you either want to download a big Study (10 000 Samples(columns) 55000 rows)
#Excel won´t be capable of loading the dataset due to its oversize
# So you have to change it here... There are 2 ways to do that

#First way: In case you just want some columns you can extract them in a basic way!


library(readr)
rcounts <- read.csv("(GTEX_data1.csv")


Just_some_columns <- rcounts[ c( "SRR1069188", "SRR1071289",
                                 "SRR1071880",
                                 "SRR1072178",
                                 "SRR1072367",
                                 "SRR1072504",
                                 "SRR1081741",
                                 "SRR1084842",
                                 "SRR1085495",
                                 "SRR1310136",
                                 "SRR1313642",
                                 "SRR1316254",
                                 "SRR1320071",
                                 "SRR1335400")]
#You will see that the Ensemble ID´s are missing, we have to fix that
Ensemblenames <- rcounts [1]



library(dplyr)
library(tidyr)
#We want the Ens_ID without the version number
Ens_ID_versionfree <-pivot_longer(Ensemblenames, cols = everything()) %>%
  separate(value, into = c('gene_id', 'beta'), sep = "\\.") %>%
  select(-name)

Ens_ID_only <- Ens_ID_versionfree[,1]

#Now we want them as our row names and so doublenames in case we had 2 versions

Ens_ID_data <- as.data.frame(Ens_ID_only)
Ens_ID_unique <-make.names(Ens_ID_data[,1], unique = TRUE)
row.names(Just_some_columns) <- Ens_ID_unique

#Thats basically it so feel free to use this for large data

#2nd way to work with huge datasets is a way that avoids excel
#It can also be used for smaller datasets but somehow doesn´t work allways
#This is why I put the version with the excel conversion on top
#It has allways worked so far

#But if the data from recount is compatible it is a faster way :-)


library(readr)
rcounts <- read.csv("(GTEX_data1.csv")

Entire_GTEX <- rcounts[-1]

#Now we do the same steps as for the first steps
Ensemblenames <- rcounts [1]



library(dplyr)
library(tidyr)
#We want the Ens_ID without the version number
Ens_ID_versionfree <-pivot_longer(Ensemblenames, cols = everything()) %>%
  separate(value, into = c('gene_id', 'beta'), sep = "\\.") %>%
  select(-name)

Ens_ID_only <- Ens_ID_versionfree[,1]

#Now we want them as our row names and so doublenames in case we had 2 versions

Ens_ID_data <- as.data.frame(Ens_ID_only)
Ens_ID_unique <-make.names(Ens_ID_data[,1], unique = TRUE)
row.names(Entire_GTEX) <- Ens_ID_unique

#Depending on the Power of your PC I recommend head() to see if it worked

head(Entire_GTEX[, c(1,2,3)], 3)

#Voila , It should be done so far, but beware really huge datasets
#like the GTEX might no be able to process in full for some PC´s 
# and packages. So in case use the first way for your interesting samples


#GTEX Annotation can be found https://gtexportal.org/home/datasets
#The GTEX recount version is V6
#look for: A de-identified, open access version of the sample annotations available in dbGaP.	GTEx_Data_V6_Annotations_SampleAttributesDS.txt
#Download it by Strg S as a txt file

# Open a tab-delimited text file with read.table()
my_file <- read.table("path/to/myfile.txt", sep="\t", header=TRUE)

# Print the contents of the file
print(my_file)



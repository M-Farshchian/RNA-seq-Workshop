############################
#####    R Basics      #####
############################
print("Hello")
print("R is fun!")

#Assigning Operator
x <- 3
y <- 5
#Variables
z <- 7
z
number <- z + y
### Tips on variable names

#-----------------
# Basic operations
#-----------------
x+y
x-y
x*y
x^y
x/y
log10(x)

### simple statistics
max(1,2,5)
min(1,2,3)
median()
mean()
range()

###help in R
help(log)
example(log)
? log

######## data type ######## 
#numeric
mode(55)
mode(4.5)
mode(-9)
#character
mode("salam")
mode("R coder")
mode('R coder')
#logical
mode(TRUE)

9 > 8
10 >16
x <- 7
x > 6
x<- 80
80 -> x
x > 6
z=19
z==18
z==19


######### Data Structures ########
## vectors
my_number<- c(1,2,3,4,5,6,7,8)
my_names<- c("cat","dog","fish")
my_logic <- c(TRUE,FALSE,TRUE)
mixed<- c(1,2,3,"ali","hossein")
new_number<- c(9,10,11)
new_vector<- c(my_number,my_names)
p <- c(1,2,3,4)
numbers1<- c(5,6,7,8)
c(p,numbers1)
v<- c(1,2,3,"ali")

# numeric vector = collection of numbers
z<-c(1,2,5,70,1000)

## we can do arithmetic operation with numeric vectors:
z*2
z/10
y<-c(4,5,6,1,9,7)
sort(y)
sum(y)
#Recycling Rule
length(y)
length(z)
# 2 vector with the same length: element_by_element computation
## what happen if we have 2 vectors with different length???
z+y

##How to make a sequence of numbers from 1 to 10 ?
b<-c(1,2,3,4,5,6,7,8,9,10)

#n:m
#n,n+1,n+2,....,m
1:10

100:80
54:32
seq(from=1,to = 100)
seq(5,100)
seq(10)
seq(10,3,-0.5)

##even numbers##
seq(from=1,to=10,by=3)
rep("Hi",time=20)

## simple plot
x <- c(21, 62, 10, 53)
names <- c("Mashhad", "Tehran", "Shiraz", "Isfahan")
pie(x,names)

install.packages("plotrix")
library(plotrix)
pie3D(x,labels = names,explode = 0.1)


####################################
##### worksapce and files #######
####################################

## working directory= where R session interacts with hard drive
## we can read data from working directory
## we can save data (plot,figure,table,..) in working directory

## absolute vs relative path
## path must be in ""
## windows vs mac path
getwd()
setwd("Desktop/workshop/data")
getwd()
list.files()



# charcter vector
my_names<-c("Ross","Rachel","Monica","Chandler")

## Matrix and data frames are rectangular data types
# They store tabular data with rows and columns
# Matrix
my_vector<- 1:20
dim(my_vector)  ##!!!
length(my_vector)

m<- matrix(data=5,nrow=2,ncol=2)
w<- matrix(my_vector,nrow=10,ncol=2)
u<- matrix(my_vector,nrow=2,ncol=10,byrow = TRUE)


# another way for making a matrix
?matrix
my_matrix2<- matrix(1:20,nrow=4,ncol=5)

##Indexing=Access elements of vector 

my_vector<-c(1,3,5,7,8,9,10,50,69,87)

# Numeric index for accessing vector elements
my_vector[3]
my_vector[c(3,5)]
my_vector[c(1,3,4,5)]
my_vector[7:10]



eng<- letters

# First element
x<-eng[1]

# Third and fourth element
eng[c(3, 4)]

# last element
eng[length(eng)] 
eng[26]
length(eng)

# Even letters
eng[seq(2, 26, 2)]

# Odd letters
eng[seq(1,25,2)]
eng[-seq(2, 26, 2)] #equivalent

## indexing by excluding elements
eng[-1]
eng[-26]
eng[-c(1,2,3,4,5)]


#### matrix= a data structure for storing objects of the same type

data <- 1:6

# Creating the matrix
matrix(data)

#you can set the number of columns or the number of rows with the ncol and nrow arguments, respectively
#you can specify if the matrix is ordered by rows or by columns with the byrow argument
# By columns
matrix(data, ncol = 2, byrow = FALSE) # byrow = FALSE by default
matrix(data, ncol = 2, nrow = 3) # Equivalent
matrix(data, nrow = 3) # Equivalent

# By rows
matrix(data, ncol = 2, byrow = TRUE)

## We can bind data tables together with cbind and rbind function
#cbind= column binding 
#rbind= row binding 
# Note that the output will be of class matrix.
x <- c(2, 7, 3, 6, 1)
y <- c(3, 7, 3, 5, 9)

# By columns
new_table1<-cbind(x, y)

# By rows
new_table2<-rbind(x, y)

# Output class
class(new_table1)
class(new_table2)

# you can use any data type inside a matrix, as long as they are homogeneous.

matrix(c(TRUE, TRUE, FALSE, TRUE), ncol = 2)
matrix(c("red", "green", "orange", "black"), ncol = 2)

#let's convert it to a table
cbind(1:6,1:4)
cbind(names,my_matrix)

## matrix see all the value as charcter
## What if??
my_names<-c("Ross","Rachel","Monica","Chandler")
my_matrix<-matrix(1:20,4,5)
cbind(my_names,my_matrix)

## How to maitain the integrity?
my_data<- data.frame(my_names,my_matrix)
class(my_data)
my_data

# add column name to a dataframe
cnames<- c("patients","age","weight","height","rate","test")
colnames(my_data)<- cnames
my_data

# add row names to a dataframe
rnames<- c("p1","p2","p3","p4")
rownames(my_data)<- rnames

##DATA FRAMES: powerful and flexible + mimic tabular dataset

##initaiting a df from row data
r1<-data.frame(a=1,b=2,c="X")
r2<-data.frame(a=3,b=4,c="Y")
r3<-data.frame(a=5,b=6,c="Z")

r<-rbind(r1,r2,r3)
class(r)
##cbind vs rbind

### Input & Output ###

#Entering data from keyboard
scores<-c(60,75,100, 94)

#Reading from a CSV FILE
library(tidyverse)
getwd()
tbl<- read.csv(file = "lncRNA.csv",header = TRUE,sep = "\t")
head(tbl)
write.csv(x,file="gene.csv")

#Subset dataframe by row and column position
df[1,3]
df[1:2,2:3]
df[-1,-c(2,3)]
df[,2]
df[3,] 
df[3,c(1,2,3)] # equivalent

## factor= represent categorical data

# we can convert character input into a factor
days <- c("Friday", "Tuesday", "Thursday", "Monday", "Wednesday", "Monday",
          "Wednesday", "Monday", "Monday", "Wednesday", "Sunday", "Saturday")
my_factor <- factor(days)
my_factor
#By default, converting a character vector to factor will order the levels alphabetically.

# we can convert numeric input into a factor
city <- c(3, 2, 1, 4, 3, 2)
my_factor2 <- factor(city)
my_factor2

#Relevel and reorder factor levels

order_days<-c("Monday","Tuesday","Wednesday","Thursday","Friday", "Saturday","Sunday")
factor(days,levels = order_days)

a<-c(2,4,67,1,0)
sort(a)

########################################################################
########################################################################
################### DIFFERNTIAL EXPRESSION ANALYSIS ####################
########################################################################
########################################################################

# Installing Required Packages

#install.packages("BiocManager")
#install.packages("devtools")
#install.packages("gplots")
#BiocManager::install("DESeq2")
#BiocManager::install("RColorBrewer")
#BiocManager::install("pheatmap")
#BiocManager::install("ggplot2")

# Loading packages
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("gplots")
library("ggplot2")
library("devtools")

# Set the working directory
directory <- "C:/Users/Reza/Desktop/Workshop/Oct 2021"
setwd(directory)

sampleFiles <- c("SRR2992186.counts","SRR2992187.counts","SRR2992194.counts","SRR2992195.counts","SRR6326445.counts","SRR6326446.counts","SRR6326447.counts","SRR6326448.counts","SRR6326449.counts")

sampleNames <- c("SRR2992186","SRR2992187","SRR2992194","SRR2992195","SRR6326445","SRR6326446","SRR6326447","SRR6326448","SRR6326449")

sampleCondition <- c("control","control","control","control","treated","treated","treated","treated","treated")

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)

sampleTable$condition <- factor(sampleTable$condition)

meta <- data.frame(sampleCondition, row.names = sampleNames)

# Creating DESeq2 object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design = ~condition)

ddsHTSeq$condition

# Setting the factor levels
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "control")

ddsHTSeq$condition

# Transform counts for data visualization
rld <- rlog(ddsHTSeq)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)

rld_cor <- cor(rld_mat)

# Plot heatmap (SDM)
pheatmap(rld_cor)

# Plot PCA 
plotPCA(rld, intgroup="condition")

# Run DESeq2 differential expression analysis
dds <- DESeq(ddsHTSeq)

# Results table will be generated using results() which will include:
# base mean, Log2 fold changes, p values and adjusted p values
res <- results(dds)

res

summary(res)

# Filter results by adjusted p value, fold change less than 1, & basemean greater than 50
res05 = subset(res, padj<0.05)

head(res05)

# Order results by padj value (the most significant to the least significant)
resOrdered <- res05[order(res05$padj),] 

write.csv(resOrdered, "results-DESeq2-normalized.csv")

# Plot counts
plotCounts(dds, gene="LINC00365", intgroup="condition")

# MA plot 
# Genes with padj < 0.1 are colored blue
plotMA(res, main = "RNAseq experiment")

# Volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-8,11)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Heatmap of the top 500 DEGs

select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:500]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=500)
heatmap.2(assay(rld)[select,], col=my_palette,
          scale="row", key=T, keysize=1,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Heatmap of the top 500 DEGs")

# Send normalized counts to tab delimited file for GSEA, etc.
gsea <- as.data.frame(counts(dds,normalized =TRUE))
gsea
write.csv(gsea, "GSEA_normalized_counts.csv")


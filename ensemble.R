setwd("~/Documents/phd/code")
AllFiles <- list.files() # This lists the files in the directory
str(AllFiles)
summary(AllFiles)
#pdf("ReflectancePlots.pdf")#Create a data file to save the plots in

setwd("~/Documents/phd/code")
whiteref <- read.table('White Standard.txt')
setwd("~/Documents/phd/code/lab treatment measurement")
AllData <- data.frame() #set up data frame of data
str(AllData) #structure of data
DataNames <- list() #list the data
str(DataNames)
topFolders <- list.files() #list the files
str(topFolders)
catagories <- data.frame()
str(catagories)
# Cycle through Folders
for(iFo in 1:length(topFolders)) {
  # Cycle through files in each folder
  folderFiles <- list.files(path = paste('./', topFolders[iFo], sep = ''))
  for(iFi in 1:length(folderFiles)) {
    dat <- read.table(paste('./', topFolders[iFo], '/', folderFiles[iFi], sep = ''), sep = '\t', header = F, skip = 15)
    dat <- apply (dat, 2, function(x) x*whiteref[,2]) #apply white reference
    AllData <- rbind(AllData, t(dat[,2]))
    DataNames <- rbind(DataNames, paste('Cat_', iFo, '_sample_', iFi, sep = ''))
    catagories <- rbind(catagories, iFo)
  }
}

rownames(AllData) <- DataNames

fit <- prcomp(AllData, scale = T) #pca

plot(fit, type = "barplot") #plot pca as bar

biplot(fit)

# Extracts for PC1 and PC2
print(fit)
pcs <- fit$x[,1:2]

plot(pcs, col = unlist(catagories))

pcd <- fit$x [,1:3] #PC1 and 3

plot (pcd, col = unlist(catagories))

pct <- fit$x[,2:3] #PC 2 and 3

plot (pct, col = unlist(catagories))

#legend("bottomright", c("P_I_M", "NMg_I_M", "K_I_M", "NKMg_I_M", "PK_I_M", "NPMg_I_M", "All_I_M", "Mg_I_M", "N_I_M", "None_I_M", "NPK_I_M", "PKMg_I_M", "All_M", "All_Control", "All_I"))

# Remove a two rows:

# fit <- prcomp(AllData[c(-3, -73)], .scale = T)

# pcs <- fit$x[,1:2]

# plot(pcs, col = unlist(catagories[c(-3, -73)]))

#
#
#
# PLS analysis
#
#
#


AllData <- data.frame()

DataNames <- list()

topFolders <- list.files()

catagories <- data.frame()
# Cycle through Folders
for(iFo in 1:length(topFolders)) {
  # Cycle through files in each folder
  folderFiles <- list.files(path = paste('./', topFolders[iFo], sep = ""))
  for(iFi in 1:length(folderFiles)) {
    dat <- read.table(paste('./', topFolders[iFo], '/', folderFiles[iFi], sep = ""), sep = "\t", header = F, skip = 15)
    AllData <- rbind(AllData, t(dat[,2]))
    DataNames <- rbind(DataNames, paste('Cat_', iFo, '_sample_', iFi, sep = ""))
    catagories <- rbind(catagories, iFo)
  }
}

rownames (AllData) <- DataNames
colnames (AllData) <- dat[,1]
em_range <- dat [,1]
setwd("~/Documents/phd/code")
source ('rse09382-mmc2.r')

pdf("PLS_Model_Fit.pdf")
setwd("~/Documents/phd/code/lab treatment measurement")
ensfit <- ensemble (AllData, unlist(catagories), em_range)
9
plot(ensfit)

stop()


scaledalldata <- scale(AllData)
library(mclust)
clustered <- Mclust(AllData)
summary(clustered)
pdf('cluster.pdf')
plot(clustered)
0
dev.off()
plot(clustered, what= 'classification')
par(mfrow = c(2,2))
plot(clustered, what = "uncertainty", dimens = c(2,1), main = "")
plot(clustered, what = "uncertainty", dimens = c(3,1), main = "")
plot(clustered, what = "uncertainty", dimens = c(2,3), main = "")
par(mfrow = c(1,1))

par(mfrow = c(2,2))
plot(clustered, what = "classification", dimens = c(2,1), main = "")
plot(clustered, what = "classification", dimens = c(3,1), main = "")
plot(clustered, what = "classification", dimens = c(2,3), main = "")
par(mfrow = c(1,1))

par(mfrow = c(2,2))
plot(clustered, what = "density", dimens = c(2,1), main = "")
plot(clustered, what = "density", dimens = c(3,1), main = "")
plot(clustered, what = "density", dimens = c(2,3), main = "")
par(mfrow = c(1,1))

ICLAllData <- mclustICL(AllData)
summary(ICLAllData)
plot(ICLAllData)
LRT = mclustBootstrapLRT(AllData, modelName = "VEI")
LRT

plot(clustered, modelName = 'VEI')
1
0
framealldata<-as.data.frame(AllData)
library(ggplot2)
View(pcs)
ndf <- data.frame(CNames = rownames(pcs), PC1 = pcs[,1], PC2 = pcs[,2])
ggplot(ndf, aes(x= PC1, y= PC2, colour="green", label= CNames))+ geom_point() +geom_text(aes(label= CNames),hjust=0, vjust=0)
head(ndf)
View(catagories)
ndf2<-ndf
ndf2[,1]<-catagories[,1]
head(ndf2)
uncat <- unlist(catagories)
ggplot(ndf2, aes(x= PC1, y= PC2, colour= 'red' , label= CNames))+ geom_point() +geom_text(aes(label= CNames),hjust=0, vjust=0, col=unlist(catagories))
ggplot(ndf2, aes(x= PC1, y= PC2, col= uncat , label= CNames)) + geom_point()
                                                                           
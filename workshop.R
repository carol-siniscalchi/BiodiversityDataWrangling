library(caper)
library(geiger)
library(stringr)
library(phytools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(diversitree)
library(phangorn)

## set up working directory, replace path with your own path
setwd("~/Desktop/workshop/")
getwd()

###########
## Fetch tree
##########

# read tree (newick) file
tree = read.tree("workshop.tree.tre")

### check if tree is ultrametric
if(is.ultrametric(tree)) {
  print("TRUE")
} else {
  tree <- force.ultrametric(tree, method="extend")
}

## plot tree
pdf(file = "tree.fullnames.pdf", width = 40, height = 160)
plotTree(tree)
dev.off()

#### QUESTION: what do the tip labels look like in this tree?

head(tree$tip.label)

#### QUESTION: should we fix the names in the tree to match environmental data?
### cleaning up tip labels on the tree
tree$tip.label <- gsub("_FMN.*", "", tree$tip.label)
tree$tip.label <- gsub("Mimosoid_", "", tree$tip.label)
tree$tip.label <- gsub("Cercidoideae_", "", tree$tip.label)
tree$tip.label <- gsub("Detarioideae_", "", tree$tip.label)
tree$tip.label <- gsub("Caesalpinioideae_", "", tree$tip.label)
tree$tip.label <- gsub("_KIB.*", "", tree$tip.label)

### remember that if we need this file in the future, we need to save it, because it's just an object right now
write.tree(tree, file = "new.workshop.tre")

########
## Diversification rates
#######

# with this dataset, we can also investigate if different environmental conditions are related to shifts in diversification.
# we'll calculate the tip diversification rates for the species on the tree

# Function from Jetz et al. (2012), aslo see Harvey et al. (2016)
DR_statistic <- function(tree, return.mean = FALSE){
  rootnode <- length(tree$tip.label) + 1
  sprates <- numeric(length(tree$tip.label))
  for (i in 1:length(sprates)){
    print(i)
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- tree$edge.length[tree$edge[,2] == node]
      node <- tree$edge[,1][tree$edge[,2] == node]			
      qx <- qx + el* (1 / 2^(index-1))			
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  if (return.mean) {
    return(mean(sprates))		
  } else {
    names(sprates) <- tree$tip.label
    return(sprates)
  }
}

#### calculate DR
DR <- DR_statistic(tree)

### save the DR in a file:
write.csv(paste0(names(DR), ",", DR, sep=""), "tip_DR_mimosoid.csv", row.names = FALSE, quote=FALSE)

### use this to plot the DR as a continuous trait in the tree
obj = contMap(tree,DR,plot=FALSE)
pdf(file="mimosoid.DR.pdf", width = 70, height = 200)
plot((setMap(obj,invert=TRUE)),legend=0.7*max(nodeHeights(tree)),sig=2,fsize=c(0.7,0.9),leg.txt="DR")
dev.off()

###########
## Fetch environmental datasets
##########

#read csv file 
BIO1 <- read.csv("BIOCLIM_1.average.csv", header = FALSE, 
                 stringsAsFactors = FALSE, row.names = NULL, sep = ",")
#give the columns names
colnames(BIO1) <- c("species", "BIO1")
#keep only distinct values?
BIO1 <- distinct(BIO1, species, .keep_all= TRUE)

#read csv file 
elevation <- read.csv("GTOPO30_ELEVATION.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
#give the columns names
colnames(elevation) <- c("species", "elevation")
#keep only distinct values?
elevation <- distinct(elevation, species, .keep_all= TRUE)

# we'll add DR as a variable in this data frame
#read csv file 
DR <- read.csv("tip_DR_mimosoid.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
#give the columns names
colnames(DR) <- c("species", "DR")
#keep only distinct values?
DR <- distinct(DR, species, .keep_all= TRUE)

## merge environmental datasets into a single dataset
combined = merge(BIO1, elevation, by = "species")
combined = merge(combined, DR, by = "species")

#check type of object
is.data.frame(combined)

#### QUESTION: what do the species labels look like in this dataframe?
head(combined)


#########
### Fetch nodulation dataset
########

#### read csv
nodulation <- read.csv("mimosoid.nod.status.csv")

##### QUESTION: what does this dataset look like? 
# How different it is from the environmental datasets? 
# How could we fix it?
head(nodulation)

# let's create a new column on our dataframe with genus names
combined$genus <- combined$species
combined$genus <- str_replace(combined$genus, "_.*", "")

# grabbing all the genus names and nodulation status
change_list <- setNames(as.character(nodulation$Nodulating), 
                        as.character(nodulation$Genus))

### apply the generic nodulation status to all species, based on the new "Genus" column 
combined$nodulation <- lapply(combined$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)

### let's get rid of rows that have unknown or NA values
combined.fixed <- na.omit(combined)
combined.fixed <- combined.fixed[!(combined.fixed$nodulation == "Unknown"),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulation == "variable"),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulation == "NULL"),]
combined.fixed <- combined.fixed[!is.null(combined.fixed$nodulation),]

### QUESTION: how many observations did we lose in the clean-up?

# let's make sure our object has the right format
combined.fixed <- as.data.frame(lapply(combined.fixed, unlist))
combined.fixed$nodulation <- as.factor(combined.fixed$nodulation)

# let's give the row names the same names as the species 
rownames(combined.fixed) <- combined.fixed$species

# let's remove the "genus" column
combined.fixed$genus <- NULL

names(combined.fixed)[names(combined.fixed) == 'x'] <- 'combined'

### ???
row.names(combined.fixed) <- combined.fixed$species

#### now let's make sure all our tree tips are matching with environmental values
#### this is a great function from the caper package
treedata_object = treedata(tree, combined.fixed)

### how many observations did we lose?

### creating new tree object and dataframe with the trimmed observations
tree.reduced <- treedata_object$phy
data.reduced <- treedata_object$data

## making sure everything is the right object type
data.reduced <- data.frame(matrix(unlist(data.reduced), nrow=nrow(data.reduced)), stringsAsFactors=FALSE)
colnames(data.reduced) = c(colnames(combined.fixed))
data.reduced$nodulation <- as.factor(data.reduced$nodulation)

##########
## Plotting some fun stuff
#########

# One question that we can ask with this dataset is if nodulators in the Mimosoid clade inhabit different
# environmental conditions than non-nodulators.
# We can use several tools for that, but first, let's plot a simple violin plot to see the differences between them.


pdf(file="BIO1.vs.nodulation.pdf", width = 6, height = 6)
ggplot(combined.fixed, aes(x = nodulation, y = BIO1, fill = nodulation)) + 
  geom_violin(trim = FALSE) + ylim(min(0), (max(300)))
dev.off()

pdf(file="elevation.vs.nodulation.pdf", width = 6, height = 6)
ggplot(combined.fixed, aes(x = nodulation, y = elevation, fill = nodulation)) + 
  geom_violin(trim = FALSE) + ylim(min(combined.fixed$elevation), (max(combined.fixed$elevation)))
dev.off()

pdf(file="DR.vs.nodulation.pdf", width = 6, height = 6)
ggplot(combined.fixed, aes(x = nodulation, y = DR, fill = nodulation)) + 
  geom_violin(trim = FALSE) + ylim(min(combined.fixed$DR), (max(combined.fixed$DR)))
dev.off()

#### now let's run a linear regression and a PGLS to test the environment vs nodulation question

## linear model
BIO1.l <- (lm(formula = BIO1 ~ nodulation, data = data.reduced))

## to check the results:``
summary(BIO1.l)

#what can we infer from the results? Is there a significant relationship?

#### PGLS -  we'll use the pgls function from the package caper, 
## for this, we need to assemble a covariance matrix of the characters in relation to their phylogenetic relationships
## this step might take a minute. When working with very large datasets, it can take a few hours, 
# so remember to save the resulting object

data_object <- comparative.data(tree.reduced, data.reduced, species, vcv=TRUE)
save(data_object, file = "combined_pgls_object.Rdata")

## now the pgls function
BIO1.p <- pgls(formula = as.numeric(BIO1) ~ nodulation, data = data_object)

## check the results:
summary(BIO1.p)

#what can we infer from the results? How is it different from the linear regression? What can we infer from this?


## let's repeat the linear regression and pgls but with DR instead of nodulation
BIO1.l.dr <- (lm(formula = BIO1 ~ DR, data = combined.fixed))
summary(BIO1.l.dr)

BIO1.p.dr <- pgls(formula = as.numeric(BIO1) ~ as.numeric(DR), data = data_object)
summary(BIO1.p.dr)



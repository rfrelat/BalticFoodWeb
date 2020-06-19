# Qualitative and quantitative description of the food web in Gulf of Riga
# Introduction to marine food web analyses
# By R. Frelat and S. Kortsch
# Last update: 22nd June 2020

# For more information, see https://rfrelat.github.io/BalticFoodWeb.html

# 1. Load and visualize the food web ------------

library(igraph)
library(fluxweb)

# Load the dataset
load("BalticFW.Rdata")

vcount(net)
ecount(net)

dim(info)
names(info)

# Visualization of the food web
# See the functional group categories 
levels(info$fg)
# Color the functional groups
colFG<- c("orange", "darkgrey", "blue", "green", "cyan")
# Assign the colors to each node (taxon)
info$colfg <- colFG[as.numeric(info$fg)]

# Plot the foodweb
# Nodes are colored (col) according to functional group 
# Links (edge) size is reduced with edge.width and edge.arrow.size
plotfw(net, col=info$colfg, 
       edge.width=0.3, edge.arrow.size=0.3)

# Overview of the food web
# basal species 
V(net)$name[degree(net, mode="in")==0]

# top predators 
V(net)$name[degree(net, mode="out")==0]

netmatrix <- get.adjacency(net, sparse=F)
heatmap(netmatrix, Rowv=NA, Colv=NA, scale="none")

# 2. Topological indicators ---------------------

# Species richness
S <- vcount(net)
S

# Connectance
C <- ecount(net)/(S*(S-1))
# same results with edge_density(net, loops=F) 
C

# Generality
# Identify predator nodes, i.e. taxa with at least one prey
pred <- degree(net, mode="in")>0
# Compute mean generality of the food web, i.e. mean number of prey per predators
G <- sum(degree(net, mode="in")[pred])/sum(pred)
G

# Vulnerability
# Identify prey nodes, i.e. taxa with at least one predator
prey <- degree(net, mode="out")>0
# Compute the mean vulnerability, i.e. mean number of predators per prey
V <- sum(degree(net, mode="out")[prey])/sum(prey)
V

# Mean shortest path
# shortest path length between all pair of nodes
sp <- shortest.paths(net)
# Mean shortest path length between different species
# [upper.tri(sp)] remove the diagonal.
# Diagonal is by default set to 0 (path to itself)
# which artificially lower the mean shortest path.
ShortPath <- mean(sp[upper.tri(sp)]) 
ShortPath

# Trophic level
# Compute the trophic level for each node  
tlnodes <- trophiclevels(net)
# Calculate the average trophic level of the food web
TL <- mean(tlnodes)
TL

# Omnivory
# Link the trophic level to the interactions
webtl <- netmatrix*tlnodes
# Remove the trophic level when no interactions
webtl[webtl==0] <- NA

#Compute the standard of the trophic levels of prey 
omninodes <- apply(webtl,2,sd, na.rm=TRUE)

# Average the standard deviation over all taxa (with more than 2 preys)  
Omni <- mean(omninodes, na.rm=TRUE)
Omni

# 3. Node-weighted indicators -------------------

#Rename the variable to shorten the R code
biomass <- info$meanB

#Define the colors of the main functional guilds
colFG <- c("orange", "darkgrey", "blue", "green", "cyan")
par(mfrow=c(1,2), mar=c(7,4,1,1))
boxplot(biomass~info$fg, las=2, col=colFG,
        ylab="Biomass (g/day/km2)", xlab="")

#Calculate percentage biomass per functional group
percB <- tapply(biomass, info$fg, sum)/sum(biomass)*100
barplot(as.matrix(percB), col=colFG, ylab="%")

#Visual representation parameter
Vscale <- 25 #multiplying factor
Vmin <- 4 #minimum size of the node
#scale the size of the node to the mean biomass
nodmax <- max(biomass)
sizeB <- (biomass/nodmax)*Vscale +Vmin

#Plot the food web
plotfw(net, col=info$colfg, size=sizeB,
       edge.width=0.3, edge.arrow.size=0.3)

# Node-weighted Connectance
nwC <- sum(degree(net)*biomass)/(2*sum(biomass)*(vcount(net)-1))
nwC

# Node-weighted generality
# Identify predators, i.e. taxa with at least one prey
pred <- degree(net, mode="in")>0
# Compute weigthed in-degree average among predators
nwG <- sum((degree(net, mode="in")*biomass)[pred])/(sum(biomass[pred]))
nwG

# Node-weighted vulnerability
# Identify prey, i.e. taxa with at least one predator
prey <- degree(net, mode="out")>0
# Compute weigthed out-degree average among prey
nwV <- sum((degree(net, mode="out")*biomass)[prey])/(sum(biomass[prey]))
nwV

#Node-weighted Trophic level
# tlnodes <- trophiclevels(net)
# Compute a weighted average of trophic levels
nwTL <- sum(tlnodes*biomass)/sum(biomass)
nwTL

# 4. Fluxweb and estimating biomass fluxes ------

# Calculate the fluxes between nodes
# netmatrix <- get.adjacency(net, sparse=F)
fluxes <- fluxing(netmatrix, biomass, info$losses, 
                  info$efficiencies, ef.level="prey")
# conversion from J/sec to kJ/day 
# 1 J/sec = 86.4 kJ/day (there are 86400 sec/day)
fluxes <- fluxes*86.4

heatmap(log(fluxes+0.00001), Rowv=NA, Colv=NA, scale="none")

# Visualize fluxes
# Create a network with link weights
netLW <- graph_from_adjacency_matrix(fluxes, weighted=TRUE)

# Set visual parameters
Escale <- 15 #multiplying coefficient
Emin <- 0.1 #minimum width 

#Calculate the width of the arrows
wid <- Emin+(E(netLW)$weight/max(E(netLW)$weight)*Escale)

# Remove the border of the frame
V(netLW)$frame.color=NA

# Plot the network
plotfw(netLW, col=info$colfg,
       edge.width=wid, edge.arrow.size=0.05)

# Link weighted indicators
fluxind(fluxes)
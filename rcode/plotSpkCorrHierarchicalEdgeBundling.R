# Use hierarchical edge bundling plot for inter-areal connections use data
# output from: createPairRelationships.m --> search for "Hierarchical Edge
# Bundling plot in R"
# Use: read.csv('[thomas-jpsth/rcode]/spkCorrVals.csv')
# Use R scripts from examples to munge data and plot

# load_libraries-------------------
# Libraries
library(dplyr)
library(ggraph)
library(igraph)
library(hashmap)

spkCorr<-read.csv('spkCorrVals.csv')
# > colnames(spkCorr)
# [1] "Pair_UID"        "X_unitNum"       "Y_unitNum"       "X_area"          "Y_area"          "XY_Dist"        
# [7] "condition"       "alignedName"     "rhoRaw_50ms"     "pvalRaw_50ms"    "rhoRaw_150ms"    "pvalRaw_150ms"  
# [13] "rhoRaw_200ms"    "pvalRaw_200ms"   "monkey"          "sessNum"         "sess"            "outcome"        
# [19] "satCondition"    "epoch"           "sameAreaPair"    "sameChannelPair" "pairCount"       "signifPlusRho"  
# [25] "signifMinusRho"  "nonSignifRho"   
# > 
# convert unitNums, sessNum to character type
cols.num <- c("X_unitNum","Y_unitNum","sessNum")
spkCorr[cols.num] <- sapply(spkCorr[cols.num],as.character)
spkCorr$signifAlpha <- spkCorr$signifPlusRho + spkCorr$signifMinusRho
spkCorr$signifAlpha[spkCorr$signifAlpha == 0] <- 0.3

# the connections are in spkCorr, which we will filter for epoch, outcome, satCondition
use.cols<-c("X_unitNum","Y_unitNum","X_area","Y_area","signifAlpha","sessNum")
connect <- select(filter(spkCorr, epoch == "PostSaccade" 
                         & outcome == "Correct" 
                         & satCondition == "Fast"
                         & X_area == "SEF"
                         & Y_area != "SEF"),
                  all_of(use.cols))
temp<-connect[connect$sessNum ==5,]

# Brain areas
d1 <- data.frame(from="origin", to=c("SEF","FEF","SC"))

xUnits <- unique(temp[,c("X_area","X_unitNum")])
colnames(xUnits)<-c("from","to")
yUnits <- unique(temp[,c("Y_area","Y_unitNum")])
colnames(yUnits)<-c("from","to")
d2 <- unique(rbind(xUnits,yUnits))
#d2 <- d2[order(d2$from),]

edges <- rbind(d1,d2)
# create a vertices data.frame. One line per object of our edges, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to)) ) )
# Let's add a column with the group of each name. It will be useful later to color points
vertices$area  <-  edges$from[ match( vertices$name, edges$to ) ]

# vertices$areaColor[vertices$area=="SEF"] <- "red"
# vertices$areaColor[vertices$area=="FEF"] <- "blue"
# vertices$areaColor[vertices$area=="SC"] <- "green"



# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( edges, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!

# The connection object must refer to the ids of the leaves:
from  <-  match( temp$X_unitNum, vertices$name)
to  <-  match( temp$Y_unitNum, vertices$name)

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", tension = .5) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  theme_void()

p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=area, size=2, alpha=0.2)) +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.4,  width=0.5, tension=0.8) +
  theme_void()
 
# Use the 'value' column of the connection data frame for the color:
p +  geom_conn_bundle(data = get_con(from = from, to = to), aes(colour=c(1:135)) )


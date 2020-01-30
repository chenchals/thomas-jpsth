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

# Brain areas
d1 <- data.frame(from="origin", to=c("SEF","FEF","SC"))
xUnits <- unique(spkCorr[,c("X_area","X_unitNum")])
colnames(xUnits)<-c("from","to")
yUnits <- unique(spkCorr[,c("Y_area","Y_unitNum")])
colnames(yUnits)<-c("from","to")
d2 <- unique(rbind(xUnits,yUnits))
edges <- rbind(d1,d2)
# create a vertices data.frame. One line per object of our edges, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to)) ) )
# Let's add a column with the group of each name. It will be useful later to color points
vertices$area  <-  edges$from[ match( vertices$name, edges$to ) ]

# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( edges, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!
monkySessNum<-unique(spkCorr[,c("monkey","sessNum")])
loop.count <- c(1:dim(monkySessNum)[1])
# the connections are in spkCorr, which we will filter for epoch, outcome, satCondition
use.cols<-c("X_unitNum","Y_unitNum","X_area","Y_area","signifAlpha","sessNum","monkey")

# temp<- select(filter(spkCorr, epoch == "PostSaccade" 
#                            & outcome == "Correct" 
#                            & satCondition == "Fast"
#                      ),)

temp<- spkCorr %>% filter(epoch == "PostSaccade" & outcome == "Correct" & satCondition == "Fast")

# SEF-FEF connections
connSefFef<-temp[temp$X_area == "SEF" & temp$Y_area == "FEF",]
# The connection object must refer to the ids of the leaves:
sefFefFrom  <-  match( connSefFef$X_unitNum, vertices$name)
sefFefTo  <-  match( connSefFef$Y_unitNum, vertices$name)
# SEF-SC connections
connSefSc<-temp[temp$X_area == "SEF" & temp$Y_area == "SC",]
# The connection object must refer to the ids of the leaves:
sefScFrom  <-  match( connSefSc$X_unitNum, vertices$name)
sefScTo  <-  match( connSefSc$Y_unitNum, vertices$name)


p<-ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=area, alpha=0.2, size=2 )) +
  theme_void()
psefFef <- p + geom_conn_bundle(data = get_con(from = sefFefFrom, to = sefFefTo), colour="red", alpha=0.9, width=0.1, tension=0.2) 
psefSc <- p + geom_conn_bundle(data = get_con(from = sefScFrom, to = sefScTo), colour="green", alpha=0.9, width=0.1, tension=0.2) 

pAll <- p + 
        geom_conn_bundle(data = get_con(from = sefFefFrom, to = sefFefTo), colour="red", alpha=0.9, width=0.5, tension=0.8) +
        geom_conn_bundle(data = get_con(from = sefScFrom, to = sefScTo), colour="green", alpha=0.9, width=0.5, tension=0.8)

sig<-temp[temp$X_area == "SEF" & temp$Y_area == "FEF" & temp$signifAlpha == 1 ,]
sefFefBndlsig<-geom_conn_bundle(data = get_con(from = match( sig$X_unitNum, vertices$name), to = match( sig$Y_unitNum, vertices$name)),
                                colour="red", alpha=0.9, width=0.5, tension=0.8)
sig<-temp[temp$X_area == "SEF" & temp$Y_area == "FEF" & temp$signifAlpha == 0 ,]
sefFefBndlnotsig<-geom_conn_bundle(data = get_con(from = match( sig$X_unitNum, vertices$name), to = match( sig$Y_unitNum, vertices$name)),
                                colour="red", alpha=0.2, width=0.5, tension=0.8)

p


for (s in loop.count)
{
  monk <- monkySessNum$monkey[s]
  sessNum <- monkySessNum$sessNum[s]
  
  temp <- spkCorr[spkCorr$sessNum == sessNum & spkCorr$monkey == monk,]
  connect <- select(filter(temp, epoch == "PostSaccade" 
                           & outcome == "Correct" 
                           & satCondition == "Fast"
                           & X_area == "SEF"
                           & Y_area != "SEF"),
                    all_of(use.cols))
  

# The connection object must refer to the ids of the leaves:
from2  <-  match( connect$X_unitNum, vertices$name)
to2  <-  match( connect$Y_unitNum, vertices$name)

# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", tension = .5) + 
#   geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
#   theme_void()

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=area, alpha=0.2, size=2 )) +
  geom_conn_bundle(data = get_con(from = from2, to = to2), colour="skyblue", alpha=0.9, width=0.1, tension=0.2) +
  geom_conn_bundle(data = get_con(from = from5, to = to5), colour="red", alpha=0.9, width=0.1, tension=0.2) +
  theme_void()
 
#p
}


# Use the 'value' column of the connection data frame for the color:
#p +  geom_conn_bundle(data = get_con(from = from, to = to), aes(colour=c(1:135)) )


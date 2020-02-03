# Use hierarchical edge bundling plot for inter-areal connections use data
# output from: createPairRelationships.m --> search for "Hierarchical Edge
# Bundling plot in R"
# Use: read.csv('[thomas-jpsth/rcode]/spkCorrVals.csv')
# Use R scripts from examples to munge data and plot

################################################################################
# Load libraries ---------------------------------------------------------------
################################################################################
library(dplyr)
library(ggraph)
library(igraph)
library(hashmap)
library(gridExtra)

################################################################################
# Define functions for plotting ---------------------------------------------
# Functions need to be compiled before they can be used
################################################################################

# fx_filter()
# Pick the connections (pairs) to be lotted so we can fine the nodeIds and do the plot 
fx_filter <- function(df,filt) 
  df %>% filter( satCondition == filt.satCond & outcome == filt.outcome & epoch == filt.epoch )

# fx_plotHierarchEdgeBundle()
# Plot the SAT correlation data as vertices in a hierarchical edge bundling layout.
fx_plotHierarchEdgeBundle <- function(df,filt,plt_base)
{
  plusColor <- "tomato"
  minusColor <- "royalblue"
  
  edgeWidth <- 0.4
  edgeTension <- 0.8
  
  # SEF-FEF correlation
  # Positive Rsc 
  temp<-df[df$X_area == "SEF" & df$Y_area == "FEF" & df$rhoRaw_150ms >= 0 & df$pvalRaw_150ms <= 0.01,]
  plt_plusFefSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                  colour=plusColor, alpha=0.8, width=edgeWidth, tension=edgeTension, linetype="solid")
  # Negative Rsc
  temp<-df[df$X_area == "SEF" & df$Y_area == "FEF" & df$rhoRaw_150ms < 0 & df$pvalRaw_150ms <= 0.01,]
  plt_minusFefSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                  colour=minusColor, alpha=0.8, width=edgeWidth, tension=edgeTension, linetype="solid")
  
  # return the plot or plot it if return val is not asked
  titleStr <- paste(toupper(filt.satCond)," - ",filt.outcome," - ",filt.epoch,sep="")
  
  plt_base + plt_minusFefSig + plt_plusFefSig +
    ggtitle(titleStr,"___ p<=0.01, - - - p>0.01")
}

################################################################################
# Do_For_AllSpkCorrs -----------------------------------------------------------
################################################################################

# Load SAT Rsc data which was pre-processed in Matlab and saved as CSV file
spkCorr <- read.csv('spkCorrVals.csv')

# Column names of spkCorr
# [1] "Pair_UID"        "X_unitNum"       "Y_unitNum"       "X_area"          "Y_area"          "XY_Dist"         "condition"       "alignedName"    
# [9] "rhoRaw_50ms"     "pvalRaw_50ms"    "rhoRaw_150ms"    "pvalRaw_150ms"   "rhoRaw_200ms"    "pvalRaw_200ms"   "X_Y_visGrade"    "X_Y_visType"    
# [17] "X_Y_moveGrade"   "X_Y_errGrade"    "X_Y_rewGrade"    "X_Y_poorIso"     "Y_visGrade"      "Y_visType"       "Y_moveGrade"     "Y_errGrade"     
# [25] "Y_rewGrade"      "Y_poorIso"       "monkey"          "sessNum"         "sess"            "outcome"         "satCondition"    "epoch"          
# [33] "sameAreaPair"    "sameChannelPair" "pairCount"       "signifPlusRho"   "signifMinusRho"  "nonSignifRho"   

# convert unitNums, sessNum to character type
cols.num <- c("X_unitNum","Y_unitNum","sessNum")
spkCorr[cols.num] <- sapply(spkCorr[cols.num], as.character)
# spkCorr$signifPlusMinus <- spkCorr$signifPlusRho + spkCorr$signifMinusRho

######################################################
# FIRST level of the hierarchy :: Area :: {SEF,FEF,SC}
######################################################

# Create a dataframe with fields "from" and "to". In this case, "from" is always set to
# 'origin', whereas "to" is set to each of our three brain areas: 'SEF','FEF','SC'.
d1 <- data.frame(from = "origin", to = c("SEF","FEF","SC"))

######################################################
# SECOND level of the hierarchy :: Neuron {1,2,3, ...}
######################################################

# Get all unique X units
xUnits <- unique(spkCorr[,c("X_area","X_unitNum")])
# For each X unit, set "area" as "from" and "unitNum" as "to"
colnames(xUnits)<-c("from","to")

# Get all unique Y units
yUnits <- unique(spkCorr[,c("Y_area","Y_unitNum")])
# For each Y unit, set "area" as "from" and "unitNum" as "to"
colnames(yUnits) <- c("from","to")

d2 <- unique(rbind(xUnits,yUnits))


# Organize edges for HEB plot
edges <- rbind(d1,d2)

# create a vertices data.frame. One line per object of our edges, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to)) ) )
# Let's add a column with the group of each name. It will be useful later to color points
vertices$area  <-  edges$from[ match( vertices$name, edges$to ) ]

# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( edges, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!
plt<-ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=area, alpha=0.2, size=2 )) +
  theme_void()
# plot Fast
filt.satCond<-"Fast"
filt.outcome<-"Correct"
filt.epoch<-"PostSaccade"
filt.sess<-"All"
p1<-fx_plotHierarchEdgeBundle(fx_filter(spkCorr,filt),filt,plt)
# plot Accurate
filt.satCond<-"Accurate"
p2<-fx_plotHierarchEdgeBundle(fx_filter(spkCorr,filt),filt,plt)

grid.arrange(p1,p2,nrow=1)





# DO_Per_Session---------------------------------------------------------------

monkySessNum <- unique(spkCorr[,c("monkey","sessNum")])
loop.count <- c(1:dim(monkySessNum)[1])
plt_list <- list()
for (s in loop.count)
{
  monk <- monkySessNum$monkey[s]
  sessNum <- monkySessNum$sessNum[s]
  
  temp <- spkCorr[spkCorr$sessNum == sessNum & spkCorr$monkey == monk,]
  
 
#p
}


# Linetypes:
#lt=c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")

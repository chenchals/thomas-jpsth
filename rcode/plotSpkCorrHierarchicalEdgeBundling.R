# Use hierarchical edge bundling plot for inter-areal connections use data
# output from: createPairRelationships.m --> search for "Hierarchical Edge
# Bundling plot in R"
# Use: read.csv('[thomas-jpsth/rcode]/spkCorrVals.csv')
# Use R scripts from examples to munge data and plot

# load_libraries---------------------------------------------------------------
# Libraries
library(dplyr)
library(ggraph)
library(igraph)
library(hashmap)
library(gridExtra)

# anonymous_functions_for_plotting---------------------------------------------
# Functions need o be compiled before they can be used...:-)
# Pick the connections (pairs) to be lotted so we can fine the nodeIds and do the plot 
fx_filter <- function(df,filt) 
  df %>% filter( satCondition == filt.satCond & outcome == filt.outcome & epoch == filt.epoch )

fx_plotIt<- function(df,filt,plt_base)
{
  # warm color
  plusColorSefSc<-"darkorange"
  plusColorSefFef<-"tomato"
  # cool color
  minusColorSefSc<-"royalblue"
  minusColorSefFef<-"slateblue1"
  w<-0.4
  t<-0.8
  # SEF-FEF Plus Rho signf and non-signif
  temp<-df[df$X_area == "SEF" & df$Y_area == "FEF" & df$rhoRaw_150ms >= 0 & df$pvalRaw_150ms <= 0.01,]
  plt_plusFefSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                   colour=plusColorSefFef, alpha=0.8, width=w, tension=t, linetype="solid")
  temp<-df[df$X_area == "SEF" & df$Y_area == "FEF" & df$rhoRaw_150ms >= 0 & df$pvalRaw_150ms > 0.01,]
  plt_plusFefNotSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                      colour=plusColorSefFef, alpha=0.2, width=w, tension=t, linetype="longdash")
  # SEF-FEF Minus Rho signf and non-signif
  temp<-df[df$X_area == "SEF" & df$Y_area == "FEF" & df$rhoRaw_150ms < 0 & df$pvalRaw_150ms <= 0.01,]
  plt_minusFefSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                    colour=minusColorSefFef, alpha=0.8, width=w, tension=t, linetype="solid")
  temp<-df[df$X_area == "SEF" & df$Y_area == "FEF" & df$rhoRaw_150ms < 0 & df$pvalRaw_150ms > 0.01,]
  plt_minusFefNotSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                       colour=minusColorSefFef, alpha=0.2, width=w, tension=t, linetype="longdash")
  
  # SEF-SC Plus Rho signf and non-signif
  temp<-df[df$X_area == "SEF" & df$Y_area == "SC" & df$rhoRaw_150ms >= 0 & df$pvalRaw_150ms <= 0.01,]
  plt_plusScSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                  colour=plusColorSefSc, alpha=0.8, width=w, tension=t, linetype="solid")
  temp<-df[df$X_area == "SEF" & df$Y_area == "SC" & df$rhoRaw_150ms >= 0 & df$pvalRaw_150ms > 0.01,]
  plt_plusScNotSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                     colour=plusColorSefSc, alpha=0.2, width=w, tension=t, linetype="longdash")
  # SEF-SC Minus Rho signf and non-signif
  temp<-df[df$X_area == "SEF" & df$Y_area == "SC" & df$rhoRaw_150ms < 0 & df$pvalRaw_150ms <= 0.01,]
  plt_minusScSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                   colour=minusColorSefSc, alpha=0.8, width=w, tension=t, linetype="solid")
  temp<-df[df$X_area == "SEF" & df$Y_area == "SC" & df$rhoRaw_150ms < 0 & df$pvalRaw_150ms > 0.01,]
  plt_minusScNotSig<-geom_conn_bundle(data = get_con(from = match( temp$X_unitNum, vertices$name), to = match( temp$Y_unitNum, vertices$name)),
                                      colour=minusColorSefSc, alpha=0.2, width=w, tension=t, linetype="longdash")
  # return the plot or plot it if return val is not asked
  titleStr<-paste(toupper(filt.satCond),"-",filt.outcome,"-",filt.epoch,sep="")
  
  plt_base + plt_minusFefNotSig + plt_minusFefSig + plt_plusFefNotSig + plt_plusFefSig +
    plt_minusScNotSig + plt_minusScSig + plt_plusScNotSig + plt_plusScSig +
    ggtitle(,"___ p<=0.01, - - - p>0.01")
  
}


# Do_For_AllSpkCorrs-----------------------------------------------------------

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
spkCorr$signifPlusMinus <- spkCorr$signifPlusRho + spkCorr$signifMinusRho

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
plt<-ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=area, alpha=0.2, size=2 )) +
  theme_void()


layout(matrix(c(1,2), byrow = TRUE))
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
p1<-fx_plotIt(fx_filter(spkCorr,filt),filt,plt)
# overwrite Fast with Accurate
filt.satCond<-"Accurate"
p2<-fx_plotIt(fx_filter(spkCorr,filt),filt,plt)

grid.arrange(p1,p2,nrow=1)


# DO_Per_Session---------------------------------------------------------------

monkySessNum<-unique(spkCorr[,c("monkey","sessNum")])
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

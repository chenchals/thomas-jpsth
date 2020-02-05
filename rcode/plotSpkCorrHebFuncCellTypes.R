# Use hierarchical edge bundling plot for inter-areal connections use data
# output from: createPairRelationships.m --> search for "Hierarchical Edge
# Bundling plot in R"
# Use: read.csv('[thomas-jpsth/rcode]/spkCorrVals.csv')
# Use R scripts from examples to munge data and plot

# load_libraries---------------------------------------------------------------
# Libraries
library(tidyverse)
library(stringr)
library(igraph)
library(ggraph)
library(gridExtra)

# anonymous_functions_for_plotting---------------------------------------------
# Functions need o be compiled before they can be used...:-)
# Pick the connections (pairs) to be lotted so we can fine the nodeIds and do the plot 
fx_filter <- function(df,filt) 
{
  df <- df %>% filter( satCondition == filt.satCond & outcome == filt.outcome & epoch == filt.epoch )
  if (!filt.sess == "All"){
    df <- df %>% filter(sess == filt.sess)
  }
  return(df)
}

fx_get_HEB <- function(df,
                       X_area,
                       Y_area,
                       plusRhoFlag,
                       sigFlag,
                       verts)
{
  w <- 0.5
  t <- 0.9
  temp <- df[df$X_area == X_area & df$Y_area == Y_area, ]
  if (plusRhoFlag == TRUE & dim(temp)[1] > 0) {
    temp <- temp[temp$rhoRaw_150ms >= 0, ]
    plusMinusLineType <- "solid"
  } else{
    temp <- temp[temp$rhoRaw_150ms < 0, ]
    plusMinusLineType <- "dotted"
  }

  if (sigFlag == TRUE & dim(temp)[1] > 0) {
    temp <- temp[temp$pvalRaw_150ms <= 0.01, ]
    a <- 0.8
    
  } else{
    temp <- temp[temp$pvalRaw_150ms > 0.01, ]
    a <- 0.2
    
  }

  if (dim(temp)[1] > 0) {
    outHeb <-
      geom_conn_bundle(
          data = get_con(from = match(temp$X_unitNum, verts$name),
          to = match(temp$Y_unitNum, verts$name)),show.legend = FALSE,
          aes(colour = visMovUnitColor),
        #aes(colour = "grey15"),
        alpha = a,
        width = w,
        tension = t,
        linetype = plusMinusLineType
      )
  } else{
    outHeb <- NULL
  }
}

fx_plotIt <- function(df, filt,  plt_base, verts)
{
  # Use call to: 
  # fx_get_HEB <- function(df,X_area,Y_area,plusRhoFlag,sigFlag,verts)
  # SEF-FEF Plus Rho signif
  plt_plusFefSig <- fx_get_HEB(df, "SEF", "FEF", TRUE, TRUE, verts)
  # SEF-FEF Minus Rho signif 
  plt_minusFefSig <- fx_get_HEB(df, "SEF", "FEF", FALSE, TRUE, verts)
  # SEF-SC Plus Rho signif 
  plt_plusScSig <- fx_get_HEB(df, "SEF", "SC", TRUE, TRUE, verts)
  # SEF-SC Minus Rho signif 
  plt_minusScSig <- fx_get_HEB(df, "SEF", "SC", FALSE, TRUE, verts)

  pltOut <- list()
  # return signif plot
  titleStr <-
    paste(toupper(filt.satCond),
          filt.outcome,
          filt.epoch,
          filt.sess,
          sep = "-")
  pltOut$signif <-
    plt_base  + plt_minusFefSig  + plt_plusFefSig + plt_minusScSig + plt_plusScSig +
    scale_edge_color_identity(guide = "none") +
    ggtitle(titleStr, "___ p<=0.01")

  return(pltOut)
  
}

# Read_spkCorr_vals_and_createEdges--------------------------------------------
# read into a tibble type, suppress column type output
spkCorr<-read_csv('spkCorrVals.csv',col_types = cols())

# > colnames(spkCorr)
# [1] "Pair_UID"        "X_unitNum"       "Y_unitNum"       "X_area"          "Y_area"          "XY_Dist"         "condition"       "alignedName"     "rhoRaw_50ms"    
# [10] "pvalRaw_50ms"    "rhoRaw_150ms"    "pvalRaw_150ms"   "rhoRaw_200ms"    "pvalRaw_200ms"   "X_visGrade"      "X_visType"       "X_moveGrade"     "X_errGrade"     
# [19] "X_rewGrade"      "X_poorIso"       "Y_visGrade"      "Y_visType"       "Y_moveGrade"     "Y_errGrade"      "Y_rewGrade"      "Y_poorIso"       "monkey"         
# [28] "sessNum"         "sess"            "outcome"         "satCondition"    "epoch"           "sameAreaPair"    "sameChannelPair" "pairCount"       "signifPlusRho"  
# [37] "signifMinusRho"  "nonSignifRho"   
# > 
# remove columns
removeCols <- c("Pair_UID", "condition", "alignedName", 
                "rhoRaw_50ms", "pvalRaw_50ms","rhoRaw_200ms", 
                "pvalRaw_200ms","sameAreaPair","sameChannelPair", "pairCount")
spkCorr <- spkCorr[ , !(colnames(spkCorr) %in% removeCols)]

# convert unitNums, sessNum to character type
numericCols <- c("X_unitNum","Y_unitNum","sessNum")
spkCorr[numericCols] <- sapply(spkCorr[numericCols], as.character)
# Find all X_ and Y_ colnames to use and remove the X_ prefix
unitColNames <- names(spkCorr)[str_detect(names(spkCorr), pattern = "X_.*")];
unitColNames <- gsub("X_","", unitColNames)
# Create allUnitNodes: get X_units
xUnits <- unique(spkCorr[,paste("X_", unitColNames, sep = "")])
# Create allUnitNodes: get X_units
yUnits <- unique(spkCorr[,paste("Y_", unitColNames, sep = "")])
# Create allUnitNodes: change column names to be same for xUnits and yUnits
colnames(xUnits) <- unitColNames
colnames(yUnits) <- unitColNames
# Create allUnitNodes: append xUnits and yUnits into a single table and get unique units 
allUnitNodes <- unique(rbind(xUnits,yUnits))
# Create allUnitNodes: Add "from" area "to" unitNum columns and polulate  
allUnitNodes[,c("from","to")] <- c(allUnitNodes[,c("area",'unitNum')])
# Create all nodes: Brain areas
brainAreas <- as_tibble(data.frame(from="brain", to=c("SEF","FEF","SC")))
# Create all nodes: add other columns
brainAreas[,unitColNames] <- NA
# Create edges
allEdges <- rbind(brainAreas,allUnitNodes)
# Create_functionalType_and_custom_sort_order----------------------------------
# Add column for vis,move,vismov types: sorting, node-color
allEdges$visMovType <- "Other"
allEdges$visMovType[abs(allEdges$visGrade) >= 2] <- "Vis"
allEdges$visMovType[abs(allEdges$moveGrade) >= 2] <- "Mov"
allEdges$visMovType[abs(allEdges$visGrade) >= 2 & abs(allEdges$moveGrade) >= 2] <- "VisMov"
# Add column for errorGarde, rewardGrade: sorting, node-shape
allEdges$errorRewardType <- "NA"
allEdges$errorRewardType[abs(allEdges$errGrade) >= 2] <- "Error"
allEdges$errorRewardType[abs(allEdges$rewGrade) >= 2] <- "Reward"
allEdges$errorRewardType[abs(allEdges$errGrade) >= 2 & abs(allEdges$rewGrade) >= 2] <- "Error & Reward"
# Add column for custom sorting and colors for visMovType column
allEdges$sortVisMov <- 0
allEdges$sortVisMov[allEdges$visMovType == "Vis"] <- 1
allEdges$sortVisMov[allEdges$visMovType == "Mov"] <- 2
allEdges$sortVisMov[allEdges$visMovType == "VisMov"] <- 3
allEdges$sortVisMov[allEdges$visMovType == "Other"] <- 4
# <chr>   <chr>      
# 1 #00BFC4 Vis        
# 2 #F8766D Mov        
# 3 #C77CFF VisMov     
# 4 #7CAE00 Other
allEdges$visMovUnitColor <- "grey50"
allEdges$visMovUnitColor[allEdges$visMovType == "Vis"] <- "#00BFC4"
allEdges$visMovUnitColor[allEdges$visMovType == "Mov"] <- "#F8766D"
allEdges$visMovUnitColor[allEdges$visMovType == "VisMov"] <- "#C77CFF"
allEdges$visMovUnitColor[allEdges$visMovType == "Other"] <- "#7CAE00"

# Add column for custom sorting for errorRewardType column
# see: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
allEdges$errorRewardTxt <- ""
allEdges$errorRewardTxt[allEdges$errorRewardType == "Error"] <- "c" #  choice error
allEdges$errorRewardTxt[allEdges$errorRewardType == "Reward"] <- "t" #  timing error
allEdges$errorRewardTxt[allEdges$errorRewardType == "Error & Reward"] <- "x" # both
# Add column for custom sorting for area column (from: column brain, SEF,FEF,SC)
allEdges$sortArea <- NA
allEdges$sortArea[allEdges$area == "SEF"] <- 0 # open square
allEdges$sortArea[allEdges$area == "FEF"] <- 1 # open circle
allEdges$sortArea[allEdges$area == "SC"] <- 2 # open triangle
# Sort allEdges
sortedEdges <- allEdges %>% group_by(sortArea,sortVisMov,errorRewardTxt)
#sortedEdges <- sortedEdges %>% arrange(sortArea,sortVisMov,visType,visGrade,moveGrade,sortErrorReward,as.numeric(unitNum),errGrade,rewGrade)
sortedEdges <- sortedEdges %>% arrange(sortArea,sortVisMov,errorRewardTxt,as.numeric(unitNum))

# Create_vertices_of_nodes_and_plot--------------------------------------------
# create a vertices data.frame. One line per object of our edges, giving features of nodes.
vertices <- as_tibble(data.frame(name = unique(c(as.character(sortedEdges$from), as.character(sortedEdges$to)) ) ))
# Let's add a column with the group of each name. It will be useful later to color points
idx <- match( vertices$name, sortedEdges$to )
vertices$area  <-  sortedEdges$area[idx]
vertices$visMovType  <-  sortedEdges$visMovType[idx]
vertices$visMovUnitColor  <-  sortedEdges$visMovUnitColor[idx]
vertices$errorRewardType  <-  sortedEdges$errorRewardType[ idx ]
vertices$errorRewardTxt  <-  sortedEdges$errorRewardTxt[ idx ]

# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( sortedEdges, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!
plt<-ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  scale_shape_discrete(solid = FALSE) +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, shape=area,
                      colour=visMovUnitColor), stroke = 1) +
  scale_color_identity(
    "Func Type",
    labels = vertices$visMovType,
    breaks = vertices$visMovUnitColor,
    guide = "legend"
  ) +
  #Change fontface. Allowed values : 1(normal),# 2(bold), 3(italic), 4(bold.italic)
  geom_text(aes(x = x*1.1, y=y*1.1,label=errorRewardTxt,colour=visMovUnitColor),
            fontface=2, show.legend = FALSE) +
  coord_equal() +
  theme_void()

# plot_for_all_sessions-----------------------
outcomes <- c("Correct", "ErrorChoice", "ErrorTiming")
sessNames <- c("All",unique(spkCorr$sess))
oPath = "../dataProcessed/analysis/spkCorr/networkPlotsColor"

for (sess in sessNames)
{
  filt.sess <- sess
  for (outcome in outcomes)
  {
    filt.outcome <- outcome
    filt.epoch <- "PostSaccade"
    outFilename <- paste(filt.outcome, filt.epoch, filt.sess, sep = "_")
    outFilename <- paste(outFilename, ".pdf", sep = "")
    
    # get plots for Fast
    filt.satCond <- "Fast"
    pltFast <-
      fx_plotIt(fx_filter(spkCorr, filt), filt, plt, vertices)
    # get plots for Accurate
    filt.satCond <- "Accurate"
    pltAccu <-
      fx_plotIt(fx_filter(spkCorr, filt), filt, plt, vertices)
    
    # return a grob...
    ZZ <- arrangeGrob(
      grobs = c(pltFast, pltAccu),
      nrow = 2,
      ncol = 1
    )
    
    # save to pdf file
    ggsave(
      outFilename,
      path = oPath,
      ZZ,
      width = 8,
      height = 10,
      units = "in"
    )
  }
}
# other_info-------
# Linetypes:
#lt=c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")

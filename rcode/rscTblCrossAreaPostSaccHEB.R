
# Use igraph library to analyze network graphs
# see kateto.net/sunbelt2019 for details
# Install_required_pkgs_and_load-----------------------------------------------
pkgs.list <- c("igraph", "network", "sna", "ggraph", "visNetwork", "threejs", 
               "networkD3", "ndtv", "htmlwidgets")
new.pkgs <- pkgs.list[!(pkgs.list %in% installed.packages()[,"Package"])]
if(length(new.pkgs)) install.packages(new.pkgs)
for (pkg in pkgs.list)
{
  library(pkg,character.only = TRUE)
}
# load_data_for_analysis-------------------------------------------------------
rscTbl <- read.csv2("rscTblCrossAreaPostSaccade.csv",sep = ",")
# filter for significant rsc
rscTbl <- rscTbl[rscTbl$signifRsc > 0,]
# code error and reward grade

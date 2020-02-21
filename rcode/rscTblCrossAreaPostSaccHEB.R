
# setwd("~/Projects/NEIWork/SATCorrelations/git/thomas-jpsth/rcode")
# Use igraph library to analyze network graphs
# see kateto.net/sunbelt2019 for details
# Install_required_pkgs_and_load-----------------------------------------------
pkgs.list <-
  c(
    "igraph",
    "network",
    "sna",
    "ggraph",
    "visNetwork",
    "threejs",
    "networkD3",
    "ndtv",
    "htmlwidgets"
  )
new.pkgs <-
  pkgs.list[!(pkgs.list %in% installed.packages()[, "Package"])]
if (length(new.pkgs))
  install.packages(new.pkgs)
for (pkg in pkgs.list)
{
  library(pkg, character.only = TRUE)
}
# Create a vector of powers of 2 (for use in conversions from binary vectors to integers). only 4 bits
powers.of.two <- 2 ^ ((4 - 1):0)
# Convert an ID (integer) to a binary vector of appropriate length
fx_decToBitVec <- function(id) {
  as.integer(rev(head(intToBits(id), 4)))
}
# Convert a binary vector of appropriate length to an ID value (integer).
fx_bitVecToDec <- function(vec) {
  as.integer(vec %*% powers.of.two)
}
# apply bitwOr to a list of ints
fx_listBitwOr <-
  function(vec) {
    res <- 0
    for (e in vec) {
      res <- bitwOr(res, e)
    }
    return(res)
  }
# convert decimal number (of bitPattern) to patternTxt
fx_decToPatternTxt <- function(patternDec)
{
  patTxt <- ""
  if (bitwAnd(8, patternDec)) 
  {
    patTxt <- paste(patTxt, "fef+", sep = "")
  }
  if (bitwAnd(4, patternDec)) 
  {
    patTxt <- paste(patTxt, "fef-", sep = "")
  }
  if (bitwAnd(2, patternDec)) 
  {
    patTxt <- paste(patTxt, "sc+", sep = "")
  }
  if (bitwAnd(1, patternDec)) 
  {
    patTxt <- paste(patTxt, "sc-", sep = "")
  }
  return(patTxt)
}

# convert decimal number (of bitPattern) to bitPattern
fx_decToPatternBits(patternDec)
{
  return(fx_decToBitVec(patternDec))
}
# load_data_for_analysis-------------------------------------------------------
rscTbl <- read.csv2("rscTblCrossAreaPostSaccade.csv", sep = ",")
# filter for significant rsc
rscTbl <- rscTbl[rscTbl$signifRsc > 0, ]
# convert unitNums to character
rscTbl$X_unitNum <- sprintf("%03d", rscTbl$X_unitNum)
rscTbl$Y_unitNum <- sprintf("%03d", rscTbl$Y_unitNum)
# add patternTxt for fef+,fef-,sc+,sc- significant rows (pair,outcome, epoch, satCond,)
rscTbl$patternTxt[rscTbl$rscSign == 1 &
                    rscTbl$Y_area == 'FEF'] <- 'fef+'
rscTbl$patternDec[rscTbl$rscSign == 1 &
                    rscTbl$Y_area == 'FEF'] <- fx_bitVecToDec(c(1, 0, 0, 0))
rscTbl$patternTxt[rscTbl$rscSign == -1 &
                    rscTbl$Y_area == 'FEF'] <- 'fef-'
rscTbl$patternDec[rscTbl$rscSign == -1 &
                    rscTbl$Y_area == 'FEF'] <- fx_bitVecToDec(c(0, 1, 0, 0))
rscTbl$patternTxt[rscTbl$rscSign == 1 &
                    rscTbl$Y_area == 'SC'] <- 'sc+'
rscTbl$patternDec[rscTbl$rscSign == 1 &
                    rscTbl$Y_area == 'SC'] <- fx_bitVecToDec(c(0, 0, 1, 0))
rscTbl$patternTxt[rscTbl$rscSign == -1 &
                    rscTbl$Y_area == 'SC'] <- 'sc-'
rscTbl$patternDec[rscTbl$rscSign == -1 &
                    rscTbl$Y_area == 'SC'] <- fx_bitVecToDec(c(0, 0, 0, 1))
# prune table columns to requied for this analysis
dropCols <-
  c(
    "X_errGrade",
    "X_rewGrade",
    "X_poorIso",
    "Y_errGrade",
    "Y_rewGrade",
    "Y_poorIso",
    "signifPlusRho",
    "signifMinusRho"
  )

rscTbl <- rscTbl[, !(names(rscTbl) %in% dropCols)]
# process_units_for_outcomeEpochSAT--------------------------------------------
outDf <- data.frame(patternDec = integer(0),
                    patternBits = integer(0),
                    patternTxt = character(""),
                    unitNum = character(""),
                    visMovType = character(""),
                    fefPlus = list(c("")),
                    fefMinus = list(),
                    scPlus = list(),
                    scMinus = list()
                    )

outcomes <- c("Correct", "ErrorChoice", "ErrorTiming")
satConditions <- c("Fast", "Accurate")
epochs <- c("PostSaccade")
for (epoch in epochs)
{
  rscTmp <- rscTbl[rscTbl$epoch == epoch, ]
  for (outcome in outcomes)
  {
    rscTmp <- rscTmp[rscTmp$outcome == outcome, ]
    for (satCondition in satConditions)
    {
      sTmp <-  rscTmp[rscTmp$satCondition == satCondition,]
      unitNums <- unique(sTmp$X_unitNum)
      
      for (unitNum in unitNums)
      {
        uTmp <- sTmp[sTmp$X_unitNum == unitNum, ]
        patternDec <- fx_listBitwOr(unique(uTmp$patternDec))
        fefPlus <- uTmp[uTmp$patternTxt == "fef+",]
        fefMinus <- uTmp[uTmp$patternTxt == "fef-",]
        scPlus <- uTmp[uTmp$patternTxt == "sc+",]
        scMinus <- uTmp[uTmp$patternTxt == "sc-",]
        
        
        
      } # unitNums
    } # satConditions
  } # outcomes
} # epochs

# other---------

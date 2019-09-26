

jpsthDirs = {
    'dataProcessed/analysis/SEF-PAPER/CHOICE_ERR_PAIRS/mat'
    'dataProcessed/analysis/SEF-PAPER/TIMING_ERR_PAIRS/mat'
    };
outDirs = {
    'dataProcessed/SEF-PAPER/CHOICE_ERR_PAIRS/r_spkCounts'
    'dataProcessed/SEF-PAPER/TIMING_ERR_PAIRS/r_spkCounts'
    };

%% Files/pairs in the dirctory and pair info
dFiles = dir(fullfile(jpsthDirs{1},'JPSTH*.mat'));
dFiles = strcat({dFiles.folder}',filesep,{dFiles.name}');
cellPairInfos = table();
%% get rasters for all alignments, for each file in directory
availConds = {{'FastErrorChoice' 'AccurateErrorChoice'}
              {'FastErrorTiming' 'AccurateErrorTiming'}};
colNames = {'rasterBins','xRasters','yRasters'};
rowNames = {'Baseline','Visual','PostSaccade'};
movingWins = [ 25, 50, 100, 150, 200];
fx_mvsum = @(rasters,win) cellfun(@(x) movsum(x,win),rasters,'UniformOutput',false);

for p = 1:1 % numel(dFiles)
    
    dat = load(dFiles{p},'-regexp','(.*ErrorChoice)|cellPairInfo');
    fastEc = dat.FastErrorChoice(rowNames,colNames);
    accuEc = dat.AccurateErrorChoice(rowNames,colNames);
    
    for w = movingWins
        movWinStr = num2str(w,'%dms');
        xMat = fx_mvsum(fastEc.xRasters,w); 
        yMat = fx_mvsum(fastEc.yRasters,w); 
        fastEc.(['ySpkCount_' movWinStr]) = fx_mvsum(fastEc.yRasters,w);
        [rho,pval] = getPearsonData(xMat,yMat);
        fastEc.(['xSpkCount_' movWinStr]) = xMat; 
        fastEc.(['ySpkCount_' movWinStr]) = yMat;
        fastEc.(['rho_pval_' movWinStr]) = [rho pval];
        
        
        
        accuEc.(['xSpkCount_' movWinStr]) = fx_mvsum(accuEc.xRasters,w); 
        accuEc.(['ySpkCount_' movWinStr]) = fx_mvsum(accuEc.yRasters,w);
        
    end
    
    
    
end



function [rho_pval,rho10,rho05,rho01] = getPearsonData(xMat,yMat)
   % Get rho, pval from matlab corr function
   % use tinv to compute crtical tval for 0.1,0.05, and 0.01 pval
   % compute the expected rho vals for the corresponding exp. pVals
   % use t = r*sqrt((n-2)/(1-r^2)) for t value 
   
   [rho,pval] = cellfun(@(x,y) corr(x,y),xMat,yMat,'UniformOutput',false);
   [rho_pval] = cellfun(@(x,y) [diag(x) diag(y)],rho,pval,'UniformOutput',false);
   n = cellfun(@(x) size(x,1),xMat);
   levels = [0.1,0,05,0.01];
   tCrit = arrayfun(@(x) tinv(levels,x),n)
   rhoCrit = sqrt((tCrit.^2)./(n-2+tCrti.^2));
   rhoCrit = num2cell(rhoCrit);
   [rho10,rho05,rho01] = deal(rhoCrit{:});
 
end
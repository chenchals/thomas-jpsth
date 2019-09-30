

jpsthDirs = {
    'dataProcessed/analysis/SEF-PAPER/jpsth/CHOICE_ERR_PAIRS/mat'
    'dataProcessed/analysis/SEF-PAPER/jpsth/TIMING_ERR_PAIRS/mat'
    };
outDirs = {
    'dataProcessed/analysis/SEF-PAPER/rSpkCounts2/mat/CHOICE_ERR_PAIRS'
    'dataProcessed/analysis/SEF-PAPER/rSpkCounts2/mat/TIMING_ERR_PAIRS'
    };
for d = 1:2
    jpsthDir = jpsthDirs{d};
    outputDir = outDirs{d};
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    
    %% Files/pairs in the dirctory and pair info
    dFiles = dir(fullfile(jpsthDirs{1},'JPSTH*.mat'));
    dFiles = strcat({dFiles.folder}',filesep,{dFiles.name}');
    cellPairInfos = table();
    %% get rasters for all alignments, for each file in directory
    availConds = {{'FastErrorChoice' 'AccurateErrorChoice'}
        {'FastErrorTiming' 'AccurateErrorTiming'}};
    % rowNames and colNames to output
    colNames = {'condition','alignedName','alignedEvent','alignedTimeWin',...
        'rasterBins','xRasters','yRasters'};
    movingWins = [ 50, 100, 200, 400];
    staticWins.Baseline = [-500 -100];
    staticWins.Visual = [50 200];
    staticWins.PostSaccade = [100 300];
    fx_mvsum = @(rasters,win) cellfun(@(x) movsum(x,win,2),rasters,'UniformOutput',false);
    
    for p = 1:1 %numel(dFiles)
        out = struct();
        cellPairInfo = load(dFiles{p},'cellPairInfo');
        cellPairInfo = cellPairInfo.cellPairInfo;
        datStruct = load(dFiles{p},'-regexp','.*Error*');
        fns = fieldnames(datStruct);
        dat = table();
        for ii = 1:numel(fns)
            fn = fns{ii};
            rowNames = datStruct.(fn).Properties.RowNames;
            tempTbl =  datStruct.(fn)(rowNames,colNames);
            tempTbl.Properties.RowNames = {};
            dat = [dat;tempTbl]; %#ok<*AGROW>
        end
        
        for w = movingWins
            movWinStr = num2str(w,'%dms');
            xMat = fx_mvsum(dat.xRasters,w);
            yMat = fx_mvsum(dat.yRasters,w);
            dat.(['xSpkCount_' movWinStr]) = xMat;
            dat.(['ySpkCount_' movWinStr]) = yMat;
            [rho_pval,dat.critRho10,dat.critRho05,dat.critRho01] = getPearsonData(xMat,yMat);
            dat.(['rho_pval_' movWinStr]) = rho_pval;
        end
        
        % do static Window spike counts
        dat.rho_pval_win = repmat(struct2cell(staticWins),numel(unique(dat.condition)),1);
        dat.xSpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.xRasters,dat.rho_pval_win,'UniformOutput',false);
        dat.ySpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.yRasters,dat.rho_pval_win,'UniformOutput',false);
        dat.rho_pval_static = getPearsonData(dat.xSpkCount_win,dat.ySpkCount_win);
        
        out.cellPairInfo = cellPairInfo;
        out.spikeCorr = dat;
        oFn = fullfile(outputDir,['rscCorr_' cellPairInfo.Pair_UID{1} '.mat']);
        saveFile(oFn,out);
    end
end


function saveFile(oFn,varData)
fprintf('Saving file %s\n',oFn);
temp = varData;
temp.lastSaved=datestr(datetime());
save(oFn,'-struct','temp');
end

function [rho_pval,critRho10,critRho05,critRho01] = getPearsonData(xMat,yMat)
% Get rho, pval from matlab corr function
% use tinv to compute crtical tval for 0.1,0.05, and 0.01 pval
% compute the critical rho vals for the pVals = 0.1,0.05,0.01 to test
% for significance
% use t = r*sqrt((n-2)/(1-r^2)) for t value

[rho,pval] = cellfun(@(x,y) corr(x,y),xMat,yMat,'UniformOutput',false);
[rho_pval] = cellfun(@(x,y) [diag(x) diag(y)],rho,pval,'UniformOutput',false);
n = cellfun(@(x) size(x,1),xMat);
% these are tied to var names rho10,rho05,rho01
levels = [0.1,0.05,0.01];
tCrit = arrayfun(@(x) tinv(levels,x),n,'UniformOutput',false);
rhoCrit = arrayfun(@(x) sqrt((tCrit{x}.^2)./(n(x)-2+tCrit{x}.^2)),(1:numel(tCrit))','UniformOutput',false);
[critRho10,critRho05,critRho01] =cellfun(@(x) deal(x(1),x(2),x(3)),rhoCrit,'UniformOutput',false);
critRho10 = cell2mat(critRho10);
critRho05 = cell2mat(critRho05);
critRho01 = cell2mat(critRho01);

end
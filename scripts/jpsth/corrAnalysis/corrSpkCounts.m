

jpsthDirs = {
    'dataProcessed/analysis/SEF-PAPER/jpsth/CHOICE_ERR_PAIRS/mat'
    'dataProcessed/analysis/SEF-PAPER/jpsth/TIMING_ERR_PAIRS/mat'
    };
wavDir = 'dataProcessed/dataset/waves2';
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
        'trialNosByCondition','alignTime','xCellSpikeTimes','yCellSpikeTimes',...
        'rasterBins','xRasters','yRasters'};
    movingWins = [ 50, 100, 200, 400];
    staticWins.Baseline = [-500 -100];
    staticWins.Visual = [50 200];
    staticWins.PostSaccade = [100 300];
    fx_mvsum = @(rasters,win) cellfun(@(x) movsum(x,win,2),rasters,'UniformOutput',false);
    
    parfor p = 1:numel(dFiles)
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
            [rho_pval,dat.critRho10,dat.critRho05,dat.critRho01] = getCorrData(xMat,yMat,'Pearson');            
            dat.(['rho_pval_' movWinStr]) = rho_pval;             
        end       
        % do static Window spike counts
        dat.rho_pval_win = repmat(struct2cell(staticWins),numel(unique(dat.condition)),1);
        dat.xSpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.xRasters,dat.rho_pval_win,'UniformOutput',false);
        dat.ySpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.yRasters,dat.rho_pval_win,'UniformOutput',false);
        dat.rho_pval_static = getCorrData(dat.xSpkCount_win,dat.ySpkCount_win,'Pearson');
        
        % get waveforms
        [dat.xWaves,dat.yWaves] = getWaveforms(wavDir,cellPairInfo,dat);
        
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

function [rho_pval,critRho10,critRho05,critRho01] = getCorrData(xMat,yMat,corrMethodStr)
    % Get rho, pval from matlab corr function

    if strcmpi(corrMethodStr,'Pearson')
        corrMethod = 'Pearson';
    elseif strcmpi(corrMethodStr,'Spearman')
        corrMethod = 'Spearman';
    elseif strcmpi(corrMethodStr,'Kendall')
        corrMethod = 'Kendall';
    end
    [rho,pval] = cellfun(@(x,y) corr(x,y,'type',corrMethod),xMat,yMat,'UniformOutput',false);
    [rho_pval] = cellfun(@(x,y) [diag(x) diag(y)],rho,pval,'UniformOutput',false);
    n = cellfun(@(x) size(x,1),xMat);
    [critRho10,critRho05,critRho01] = getCriticalTvalue(n);
end

function [critRho10,critRho05,critRho01] = getCriticalTvalue(sampleSizeArray)
    % use tinv to compute crtical tval for 0.1,0.05, and 0.01 pval
    % compute the critical rho vals for the pVals = 0.1,0.05,0.01 to test
    % for significance
    % use t = r*sqrt((n-2)/(1-r^2)) for t value
    n = sampleSizeArray;
    % these are tied to var names rho10,rho05,rho01
    levels = [0.1,0.05,0.01];
    tCrit = arrayfun(@(x) tinv(levels,x),n,'UniformOutput',false);
    rhoCrit = arrayfun(@(x) sqrt((tCrit{x}.^2)./(n(x)-2+tCrit{x}.^2)),(1:numel(tCrit))','UniformOutput',false);
    [critRho10,critRho05,critRho01] =cellfun(@(x) deal(x(1),x(2),x(3)),rhoCrit,'UniformOutput',false);
    critRho10 = cell2mat(critRho10);
    critRho05 = cell2mat(critRho05);
    critRho01 = cell2mat(critRho01);
end

function [xWaves,yWaves] = getWaveforms(wavDir,cellPairInfo,dat)
    %% inline function for aligning, windowing TS of trials
    fx_alignTs = @(tsByTrl,alinTimeByTrl) arrayfun(@(idx) tsByTrl{idx} - alinTimeByTrl(idx),(1:numel(alinTimeByTrl))','UniformOutput',false);
    fx_alindTsToWin = @(alindTsByTrl, tsWin) cellfun(@(x) x(x>=tsWin(1) & x<=tsWin(2)),alindTsByTrl,'UniformOutput',false);
    % find matching indices for 2 arrrays of aligned ts
    % max abs. diff is 1 ms
    diffMax = 1; % 1 ms different
    fx_matchedSpkIdx = @(alindUTs,alindWfTs) arrayfun(@(uTs) find(abs(alindWfTs - uTs)<=diffMax,1),alindUTs,'UniformOutput',false);

    %% get all waves parse for sel. trials and match wave Ts.
    % X Unit
    unitName = sprintf('Unit_%03d',cellPairInfo.X_unitNum);
    allWavs = load(fullfile(wavDir,[unitName '.mat']));
    wavTsX = allWavs.(unitName).wavSearchTs{1};
    wavsX = allWavs.(unitName).wavSearch{1};
    % Y Unit
    unitName = sprintf('Unit_%03d',cellPairInfo.Y_unitNum);
    allWavs = load(fullfile(wavDir,[unitName '.mat']));
    wavTsY = allWavs.(unitName).wavSearchTs{1};
    wavsY = allWavs.(unitName).wavSearch{1};

    xWaves = cell(12,1);
    yWaves = cell(12,1);

    for jj = 1:size(dat,1)
        selTrls = find(dat.trialNosByCondition{jj});
        alignTime = dat.alignTime{jj};
        alignWin = dat.alignedTimeWin{jj};
        % align ts from waveform data
        alindWfTsX = fx_alignTs(wavTsX(selTrls),alignTime);
        alindWfTsXWin = fx_alindTsToWin(alindWfTsX,alignWin);
        alindWfTsY = fx_alignTs(wavTsY(selTrls),alignTime);
        alindWfTsYWin = fx_alindTsToWin(alindWfTsY,alignWin);
        % aligned windowd unit ts
        alindUnitTsXWin = dat.xCellSpikeTimes{jj};
        alindUnitTsYWin = dat.yCellSpikeTimes{jj};
        % find matching indices for extracting waveforms
        spkIdxByTrlX = cellfun(@(aUts,aWts) fx_matchedSpkIdx(aUts,aWts),alindUnitTsXWin,alindWfTsXWin,'UniformOutput',false);
        xWaves{jj} = cellfun(@(wfByTrl,idxByTrl) wfByTrl(cell2mat(idxByTrl'),:),wavsX(selTrls),spkIdxByTrlX,'UniformOutput',false);
        spkIdxByTrlY = cellfun(@(aUts,aWts) fx_matchedSpkIdx(aUts,aWts),alindUnitTsYWin,alindWfTsYWin,'UniformOutput',false);
        yWaves{jj} = cellfun(@(wfByTrl,idxByTrl) wfByTrl(cell2mat(idxByTrl'),:),wavsY(selTrls),spkIdxByTrlY,'UniformOutput',false);

    end


end


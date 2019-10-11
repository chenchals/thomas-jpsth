
%%
monkIdsToDo = {'D','E'};
jpsthDirs = {
    'dataProcessed/analysis/JPSTH-10ms/jpsth_SEF-SEF/mat' 
    'dataProcessed/analysis/JPSTH-10ms/jpsth_SEF-FEF/mat'    
    'dataProcessed/analysis/JPSTH-10ms/jpsth_SEF-SC/mat'     
    'dataProcessed/analysis/JPSTH-10ms/jpsth_FEF-FEF/mat'    
    'dataProcessed/analysis/JPSTH-10ms/jpsth_FEF-SC/mat'     
    'dataProcessed/analysis/JPSTH-10ms/jpsth_SC-SC/mat'      
    'dataProcessed/analysis/JPSTH-10ms/jpsth_SEF-NSEFN/mat'  
    'dataProcessed/analysis/JPSTH-10ms/jpsth_FEF-NSEFN/mat'  
    'dataProcessed/analysis/JPSTH-10ms/jpsth_SC-NSEFN/mat'   
    'dataProcessed/analysis/JPSTH-10ms/jpsth_NSEFN-NSEFN/mat'
    };
wavDir = 'dataProcessed/dataset/waves';
outDirs = regexprep(jpsthDirs,'JPSTH-10ms/jpsth_','spkCorr/spkCorr_');

%% do only for Da and Eu
jpsthPairsDaEu = load('dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat');
jpsthPairsDaEu = jpsthPairsDaEu.JpsthPairCellInfoDB;
jpsthPairsDaEu = jpsthPairsDaEu{ismember([jpsthPairsDaEu.X_monkey],monkIdsToDo),{'Pair_UID'}};

fx_saveFilenameOnError = @(fn) save(fn);
%%
for d = 1:numel(jpsthDirs)
    jpsthDir = jpsthDirs{d};
    outputDir = outDirs{d};
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    
    %% Files/pairs in the dirctory and pair info
    dFiles = dir(fullfile(jpsthDir,'*.mat'));
    dFiles = strcat({dFiles.folder}',filesep,{dFiles.name}');
    % restrict to Da, Eu pairs
    dFiles = dFiles(contains(dFiles,jpsthPairsDaEu));
    
    %% get rasters for all alignments, for each file in directory
    %datStruct = load(dFiles{p},'-regexp','Accurate*|Fast*');
%     availConds = {{'FastErrorChoice' 'AccurateErrorChoice'}
%         {'FastErrorTiming' 'AccurateErrorTiming'}
%         {'FastCorrect' 'AccurateCorrect'}
%         };
    %% Baseline epoch to extract from Visual epoch
    baselineWin = [-600 50];
    
    %% Static r_sc
    movingWins = [50, 100, 200, 400];
    staticWins.Baseline = [-500 -100];
    staticWins.Visual = [50 250];
    staticWins.PostSaccade = [0 400];
    staticWins.PostReward = [0 600];
    fx_mvsum = @(rasters,win) cellfun(@(x) movsum(x,win,2),rasters,'UniformOutput',false);
    fx_zscoreTrls = @(matCellArr) cellfun(@(x) zscore(x,0,2),matCellArr,'UniformOutput',false);

%    dFiles{1} = 'dataProcessed/analysis/JPSTH-5ms/jpsth_FEF-FEF/mat/JPSTH-PAIR_0697.mat';
    %% For each file
    parfor p = 1:numel(dFiles)
        out = struct();
        colNames = {'condition','alignedName','alignedEvent','alignedTimeWin',...
            'trialNosByCondition','alignTime','xCellSpikeTimes','yCellSpikeTimes',...
            'rasterBins','xRasters','yRasters','firstSortByName','firstSortByTime'...
            };
        cellPairInfo = load(dFiles{p},'cellPairInfo');
        cellPairInfo = cellPairInfo.cellPairInfo;
        %datStruct = load(dFiles{p},'-regexp','.*Error*');
        datStruct = load(dFiles{p},'-regexp','Accurate*|Fast*');
        %% process, now that we have added the baseline row (Not all cols are valid       
        fns = fieldnames(datStruct);        
        dat = table();
        for ii = 1:numel(fns)
            fn = fns{ii};
            if isempty(datStruct.(fn))
                continue
            end
            rowNames = datStruct.(fn).Properties.RowNames;            
            datFns = datStruct.(fn).Properties.VariableNames;
             tempTbl =  datStruct.(fn)(rowNames,colNames);
             tempTbl.Properties.RowNames = {};
             dat = [dat;tempTbl]; %#ok<*AGROW>
        end
        xMatRaw = dat.xRasters;
        yMatRaw = dat.yRasters;
        % Z-Score each trial
        % see: https://www.nature.com/articles/s41593-019-0477-1
        % Ruff & Cohen Simultaneous multi-area recordings suggest that
        % attention improves performance by reshaping stimulus
        % representations 2019, Nature Neuroscience 22:1669-1676
        [xMatZ,xMatMean,xMatStd] = fx_zscoreTrls(xMatRaw);
        [yMatZ,yMatMean,yMatStd] = fx_zscoreTrls(yMatRaw);
        dat.xRasters_Z = xMatZ;
        dat.xRastersMean = xMatMean;
        dat.xRastersStd = xMatStd;
        dat.yRasters_Z = yMatZ;
        dat.yRastersMean = yMatMean;
        dat.yRastersStd = yMatStd;
        
        %%
        for w = movingWins
            movWinStr = num2str(w,'%dms');
            xMat = fx_mvsum(xMatRaw,w);
            yMat = fx_mvsum(yMatRaw,w);
            dat.(['xSpkCount_' movWinStr]) = xMat;
            dat.(['ySpkCount_' movWinStr]) = yMat;
            [rho_pval,dat.critRho10,dat.critRho05,dat.critRho01] = getCorrData(xMat,yMat,'Pearson');            
            dat.(['rho_pval_' movWinStr]) = rho_pval;      
            % do on Zscored            
            xMat = fx_mvsum(xMatZ,w);
            yMat = fx_mvsum(yMatZ,w);
            dat.(['xSpkCount_' movWinStr '_Z']) = xMat;
            dat.(['ySpkCount_' movWinStr '_Z']) = yMat;
            [rho_pval,dat.critRho10_Z,dat.critRho05_Z,dat.critRho01_Z] = getCorrData(xMat,yMat,'Pearson');
            dat.(['rho_pval_' movWinStr '_Z']) = rho_pval;           
        end       
        
        %% do static Window spike counts
        dat.rho_pval_win = repmat(struct2cell(staticWins),numel(unique(dat.condition)),1);
        dat.xSpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.xRasters,dat.rho_pval_win,'UniformOutput',false);
        dat.ySpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.yRasters,dat.rho_pval_win,'UniformOutput',false);
        dat.rho_pval_static = getCorrData(dat.xSpkCount_win,dat.ySpkCount_win,'Pearson');
        
        %% do static Window spike counts - Z-scored
        dat.rho_pval_win_Z = repmat(struct2cell(staticWins),numel(unique(dat.condition)),1);
        dat.xSpkCount_win_Z = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.xRasters_Z,dat.rho_pval_win_Z,'UniformOutput',false);
        dat.ySpkCount_win_Z = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
            dat.rasterBins,dat.yRasters_Z,dat.rho_pval_win_Z,'UniformOutput',false);
        dat.rho_pval_static_Z = getCorrData(dat.xSpkCount_win_Z,dat.ySpkCount_win_Z,'Pearson');
        
        %% get waveforms
        try
        [dat.xWaves,dat.yWaves] = getWaveforms(wavDir,cellPairInfo,dat);
        dat.xWaveWidths = getWaveformWidths(dat.xWaves);
        dat.yWaveWidths = getWaveformWidths(dat.yWaves);
        catch me
            msg = sprintf('Error in processing pair_UID %s while call to getWaveforms ...\n',cellPairInfo.Pair_UID{1});
            fprintf(msg);
            getReport(me);
            oFn = fullfile(outputDir,['ERROR_PROCESSING_rscCorr_' cellPairInfo.Pair_UID{1} '.mat']);
            fx_saveWorkspaceOnError(oFn);
            continue        
        end
        
        %% save datta file for r-spkCounts
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

    xWaves = cell(size(dat,1),1);
    yWaves = cell(size(dat,1),1);

    %% get all waves parse for sel. trials and match wave Ts.
    % X Unit
    unitName = sprintf('Unit_%03d',cellPairInfo.X_unitNum);
    if ~exist(fullfile(wavDir,[unitName '.mat']),'file')
        return;
    end
    allWavs = load(fullfile(wavDir,[unitName '.mat']));
    wavTsX = allWavs.(unitName).wavSearchTs{1};
    wavsX = allWavs.(unitName).wavSearch{1};
    % Y Unit
    unitName = sprintf('Unit_%03d',cellPairInfo.Y_unitNum);
    allWavs = load(fullfile(wavDir,[unitName '.mat']));
    wavTsY = allWavs.(unitName).wavSearchTs{1};
    wavsY = allWavs.(unitName).wavSearch{1};

    for jj = 1:size(dat,1)
        selTrls = dat.trialNosByCondition{jj};
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


function [wavWidths] = getWaveformWidths(wavforms)
    % given a matrix of waveforms, find width by interpolating with a
    % cubic-spline at 10x sampling
    fx_spkWid = @(spks,magFactor) ...
        arrayfun(@(s) (find(spks(s,:)==max(spks(s,:),[],2),1) ...
                    - find(spks(s,:)==min(spks(s,:),[],2),1))/magFactor,...
                      (1:size(spks,1))');
    wavWidths = cell(numel(wavforms),1);
    interpolateAt = 10; % 10x
    for ro = 1:numel(wavforms)
        w = wavforms{ro};
        if isempty(w)
            continue;
        end      
        x = size(w{1},2);
        x = (1:x)';
        xq = (1:1/interpolateAt:max(x))';
        % there could be trials with no spikes
        notEmptyIdx = find(~cellfun(@isempty,w));
        % interpolated waveforms for trials with spikes
        wq = cellfun(@(s) interp1(x,s',xq,'pchip')',w(notEmptyIdx),'UniformOutput',false);
        % waveform widths in sampling space
        wqWideTemp = cellfun(@(s) fx_spkWid(s,interpolateAt),wq,'UniformOutput',false);
        % put back widths for correct trials
        wqWide = cell(size(w,1),1);
        wqWide(notEmptyIdx) = wqWideTemp;
        wavWidths{ro} = wqWide;
        
    end
end


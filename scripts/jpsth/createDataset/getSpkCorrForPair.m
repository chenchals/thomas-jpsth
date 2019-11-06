function [outSpkCorr] = getSpkCorrForPair(cellPair,xSpkTimes,ySpkTimes,evntTimes,trialTypes,...
    conditions,alignNames,alignEvents,alignTimeWin,firstSortEvent,wavDir)
%GETSPKCORRFORPAIR Summary of this function goes here
% cellPair : the pair row in the table JPSTH_PAIRS_CellInfoDB.mat
% xSpkTimes : X-axis cell: spike times as cell array of nTrials
% ySpkTimes : Y-axis cell: spike times as cell array of nTrials
% evntTimes : event times for session from TrialEventTimesDB.mat
% trialTypes : trial types for session from TrialTypesDB.mat
% conditions : Trial type names: AccurateCorrect, FastCorrect... from
%              TrialTypesDB
%   Example:  conditions = {
%            'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
%            'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
%             };
% alignNames : Name of the aligned event window : Baseline, Visual,
%              postSaccade
% alignEvents : Name of the event to align on ... from
%               TrialEventTimesDB.mat
% alignTimeWin : A cell array of time windows for epochs corres[ponding to
%                align events

    %% Dynamic spk corr across the trial
    movingWins = [50, 100, 200, 10];
    % Rasters are got from jpsth data:
    % alignTimeWin = {[-600 100],[-200 400],[-100 500],[-200 700]};
    % when spkCorr is computed: for each moving window:
    % Fill NaN for window times min(alignWin)-movingWin/2
    % Fill NaN for window times min(alignWin)-movingWin/2
    fx_mvsum = @(rasters,win) movsum(double(rasters),win,2,'Endpoints','fill');
    %% Static spk corr windows for computing spike corr
    staticWins.Baseline = [-500 -100];%[-500 -100];
    staticWins.Visual = [50 200];%[50 250];
    staticWins.PostSaccade = [100 300];%[0 400];
    staticWins.PostReward = [100 300];%[0 600];

    warning('off');
    % ignore processing if the sel. trials are below thisNum.
    nTrialsThreshold = 5;

    units = struct();
    spkCorr = struct();

    XCellId = ['DSP' cellPair.X_unit{1}];
    YCellId = ['DSP' cellPair.Y_unit{1}];

    units.(XCellId) = xSpkTimes;
    units.(YCellId) = ySpkTimes;
    
    units.('X_trialMeanStd') = SpikeUtils.getTrialMeanStd(xSpkTimes);
    units.('Y_trialMeanStd') = SpikeUtils.getTrialMeanStd(ySpkTimes);
    
    baselineWinForMeanStd = [-600 -100];
    temp = SpikeUtils.alignSpikeTimes(xSpkTimes,evntTimes.CueOn{1},baselineWinForMeanStd);
    temp = SpikeUtils.rasters(temp,baselineWinForMeanStd);
    temp = temp.rasters;
    units.('X_baselineMeanStd') = arrayfun(@(x) [mean(temp(x,:)) std(temp(x,:))],(1:size(temp,1))','UniformOutput',false );
    temp = SpikeUtils.alignSpikeTimes(ySpkTimes,evntTimes.CueOn{1},baselineWinForMeanStd);
    temp = SpikeUtils.rasters(temp,baselineWinForMeanStd);
    temp = temp.rasters;
    units.('Y_baselineMeanStd') = arrayfun(@(x) [mean(temp(x,:)) std(temp(x,:))],(1:size(temp,1))','UniformOutput',false );
    clear temp;
    for cond = 1:numel(conditions)
        try
            % incase something breaks continue...
            condition = conditions{cond};
            selTrials = trialTypes.(condition){:};
            if isempty(selTrials)
                spkCorr.(condition) = [];
                continue;
            end

            %% Mutually Exclusive trials for Choice/Timing Errors
            % If condition is *ChoiceErr or *TimingError ensure mutually
            % exclusive
            otherCondition = [];
            if contains(condition,'ChoiceError')
                otherCondition = regexprep(codition,'ChoiceError','TimingError');
            elseif contains(condition, 'TimingError')
                otherCondition = regexprep(codition,'TimingError','ChoiceError');
            end
            if ~isempty(otherCondition)
                selTrials(trialTypes.(otherCondition){:}) = 0;
            end

            % first check for nTrials
            if isempty(selTrials) || numel(selTrials) < nTrialsThreshold
                spkCorr.(condition) = [];
                continue;
            end

            %% check if trials need to be dropped due to poor_isolation..
            trRem = cellPair.X_trRemSAT{1};
            if ~isempty(trRem)
                selTrials(trRem(1):trRem(2)) = 0;
            end
            trRem = cellPair.Y_trRemSAT{1};
            if ~isempty(trRem)
                selTrials(trRem(1):trRem(2)) = 0;
            end
            if isempty(selTrials)
                spkCorr.(condition) = [];
                continue;
            end
            % second check for nTrials
            selTrials = find(selTrials);
            if numel(selTrials) <= nTrialsThreshold
                spkCorr.(condition) = [];
                continue;
            end

            %% Sort selected trials based on the event
            selTrlsTbl = table();
            selTrlsTbl.selTrials = selTrials;
            %% for each aligned event
            opts = struct();
            for evId = 1:numel(alignEvents)
                alignedEvent = alignEvents{evId};
                alignedTimeWin = alignTimeWin{evId};
                alignedName = alignNames{evId};
                if isempty(firstSortEvent{evId})
                    firstSortName = [];
                    selTrlsTbl.firstSortTime = -inf(size(selTrlsTbl,1),1);
                else
                    firstSortName = firstSortEvent{evId};
                    firstSortTime = double(evntTimes.(firstSortName){1});
                    temp = firstSortTime(selTrlsTbl.selTrials);
                    temp(temp==0 | isnan(temp)) = -inf;
                    selTrlsTbl.firstSortTime = temp;
                end
                selTrlsTblSorted = sortrows(selTrlsTbl,{'firstSortTime'});
                selTrialsSorted = selTrlsTblSorted.selTrials;
                alignTime = evntTimes.CueOn{1};
                if ~strcmp(alignedEvent,'CueOn')
                    alignTime = alignTime + double(evntTimes.(alignedEvent){1}(:));
                end
                %% Align Spike times and get rasters
                alignTime = alignTime(selTrialsSorted);
                XAligned = SpikeUtils.alignSpikeTimes(units.(XCellId)(selTrialsSorted),alignTime, alignedTimeWin);
                YAligned = SpikeUtils.alignSpikeTimes(units.(YCellId)(selTrialsSorted),alignTime, alignedTimeWin);
                tempRast = SpikeUtils.rasters(XAligned,alignedTimeWin);
                XRasters = tempRast.rasters;
                tempRast = SpikeUtils.rasters(YAligned,alignedTimeWin);
                YRasters = tempRast.rasters;
                rasterBins = tempRast.rasterBins;                               
                clear tempRast;
                %% Z scores
                % get Z score for given sample using means and stds specified for each trial
                % Z-Score using baseline mean and std
                fx_getZscore = @(rast,meanStd) cell2mat(arrayfun(@(x) (double(rast(x,:))- meanStd{1}(1))./meanStd{1}(2), (1:size(rast,1))','UniformOutput',false));
                xRasters_Z_baseline = fx_getZscore(XRasters,units.X_baselineMeanStd);
                yRasters_Z_baseline = fx_getZscore(YRasters,units.Y_baselineMeanStd);
                xRasters_Z_trial = fx_getZscore(XRasters,units.X_trialMeanStd);
                yRasters_Z_trial = fx_getZscore(YRasters,units.Y_trialMeanStd);

                %% gather vars fo corrSpk_PAIR_xxxx.mat
                opts(evId,1).condition = {condition};
                opts(evId,1).alignedName = {alignedName};
                opts(evId,1).alignedEvent = {alignedEvent};
                opts(evId,1).alignedTimeWin = {alignedTimeWin};
                opts(evId,1).alignTime = {alignTime};
                opts(evId,1).firstSortByName = {firstSortName};
                opts(evId,1).firstSortByTime = selTrlsTblSorted.firstSortTime;
                opts(evId,1).trialNosByCondition = selTrialsSorted;
                opts(evId,1).xCellSpikeTimes = {XAligned}; %#ok<*AGROW>
                opts(evId,1).yCellSpikeTimes = {YAligned};
                opts(evId,1).xBaselineMeanStd = {cell2mat(units.X_baselineMeanStd)};
                opts(evId,1).yBaselineMeanStd = {cell2mat(units.Y_baselineMeanStd)};
                opts(evId,1).xTrialMeanStd = {cell2mat(units.X_trialMeanStd)};
                opts(evId,1).yTrialMeanStd = {cell2mat(units.Y_trialMeanStd)};
                opts(evId,1).rasterBins = {rasterBins};
                opts(evId,1).xRasters = {XRasters};
                opts(evId,1).yRasters = {YRasters};
                opts(evId,1).xRasters_Z_baseline = {xRasters_Z_baseline};
                opts(evId,1).yRasters_Z_baseline = {yRasters_Z_baseline};
                opts(evId,1).xRasters_Z_trial = {xRasters_Z_trial};
                opts(evId,1).yRasters_Z_trial = {yRasters_Z_trial};
                
                %% Compute spike corrs across trial using different moving windows for each aligned event
                for ww = 1:numel(movingWins)
                    w = movingWins(ww);
                    movWinStr = num2str(w,'%dms');
                    % spike corrs for raw spk count
                    xMat1 = fx_mvsum(XRasters,w);
                    yMat1 = fx_mvsum(YRasters,w);
                    opts(evId,1).(['xSpkCount_' movWinStr]) = xMat1;
                    opts(evId,1).(['ySpkCount_' movWinStr]) = yMat1;
                    [rho_pval1,opts(evId,1).critRho10,opts(evId,1).critRho05,opts(evId,1).critRho01] = getCorrData(xMat1,yMat1,'Pearson');
                    opts(evId,1).(['rho_pval_' movWinStr]) = rho_pval1;
                    % spike corrs for Z scored (using Baseline mean/std) counts
                    xMat2 = fx_mvsum(xRasters_Z_baseline,w);
                    yMat2 = fx_mvsum(yRasters_Z_baseline,w);
                    opts(evId,1).(['xSpkCount_' movWinStr '_Z_baseline']) = xMat2;
                    opts(evId,1).(['ySpkCount_' movWinStr '_Z_baseline']) = yMat2;
                    [rho_pval2,opts(evId,1).critRho10_Z_baseline,opts(evId,1).critRho05_Z_baseline,opts(evId,1).critRho01_Z_baseline] = getCorrData(xMat2,yMat2,'Pearson');
                    opts(evId,1).(['rho_pval_' movWinStr '_Z_baseline']) = rho_pval2;
                    % spike corrs for Z scored (using Whole Trial mean/std) counts
                    xMat3 = fx_mvsum(xRasters_Z_trial,w);
                    yMat3 = fx_mvsum(yRasters_Z_trial,w);
                    opts(evId,1).(['xSpkCount_' movWinStr '_Z_trial']) = xMat3;
                    opts(evId,1).(['ySpkCount_' movWinStr '_Z_trial']) = yMat3;
                    [rho_pval3,opts(evId,1).critRho10_Z_trial,opts(evId,1).critRho05_Z_trial,opts(evId,1).critRho01_Z_trial] = getCorrData(xMat3,yMat3,'Pearson');
                    opts(evId,1).(['rho_pval_' movWinStr '_Z_trial']) = rho_pval3;
                end
                
                %% Compute spike corrs for static windows for each aligned event
                % Static windows spike corr for Raw counts
                opts(evId,1).rho_pval_win = {staticWins.(alignedName)};
                opts(evId,1).xSpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{XRasters},opts(evId,1).rho_pval_win,'UniformOutput',false);
                opts(evId,1).ySpkCount_win = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{YRasters},opts(evId,1).rho_pval_win,'UniformOutput',false);
                opts(evId,1).rho_pval_static = {getCorrData(opts(evId,1).xSpkCount_win{1},opts(evId,1).ySpkCount_win{1},'Pearson')};
                
                % Static windows spike corr for - Z-scored (using Baseline mean/std) count
                opts(evId,1).rho_pval_win_Z_baseline = {staticWins.(alignedName)};
                opts(evId,1).xSpkCount_win_Z_baseline = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{xRasters_Z_baseline},opts(evId,1).rho_pval_win_Z_baseline,'UniformOutput',false);
                opts(evId,1).ySpkCount_win_Z_baseline = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{yRasters_Z_baseline},opts(evId,1).rho_pval_win_Z_baseline,'UniformOutput',false);
                opts(evId,1).rho_pval_static_Z_baseline = {getCorrData(opts(evId,1).xSpkCount_win_Z_baseline{1},opts(evId,1).ySpkCount_win_Z_baseline{1},'Pearson')};

                % Static windows spike corr for - Z-scored (using Baseline mean/std) count
                opts(evId,1).rho_pval_win_Z_trial = {staticWins.(alignedName)};
                opts(evId,1).xSpkCount_win_Z_trial = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{xRasters_Z_trial},opts(evId,1).rho_pval_win_Z_baseline,'UniformOutput',false);
                opts(evId,1).ySpkCount_win_Z_trial = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{yRasters_Z_trial},opts(evId,1).rho_pval_win_Z_baseline,'UniformOutput',false);
                opts(evId,1).rho_pval_static_Z_trial = {getCorrData(opts(evId,1).xSpkCount_win_Z_trial{1},opts(evId,1).ySpkCount_win_Z_trial{1},'Pearson')};
                %% get waveforms
                [opts(evId,1).xWaves,opts(evId,1).yWaves] = getWaveforms(wavDir,cellPair,opts(evId,1));
                opts(evId,1).xWaveWidths = getWaveformWidths(opts(evId,1).xWaves);
                opts(evId,1).yWaveWidths = getWaveformWidths(opts(evId,1).yWaves);
                
                
            end % for alignEvents
            %tempSpkCorr.Properties.RowNames = alignNames;
            spkCorr.(condition) = struct2table(opts,'AsArray',true);
        catch mE
            getReport(mE)
            continue
        end
    end % for conditions
    %% Save for the current pair
    fns = fieldnames(spkCorr);
    tempTbl = table();
    for fn = 1:numel(fns)
        if ~isempty(spkCorr.(fns{fn}))
            tempTbl = [tempTbl;spkCorr.(fns{fn})];
        end
    end
    % convert all doubles to single to save space
    fns = tempTbl.Properties.VariableNames;
    for ii = 1:numel(fns)
        fn = fns{ii};
        if iscell(tempTbl.(fn)) && isa(tempTbl.(fn){1},'double')
            tempTbl.(fn) = cellfun(@(x) single(x),tempTbl.(fn),'UniformOutput',false);
        elseif isa(tempTbl.(fn)(1),'double')
            tempTbl.(fn) = single(tempTbl.(fn));
        end
    end

    outSpkCorr.spkCorr = tempTbl;
    outSpkCorr.cellPairInfo = cellPair;
    outSpkCorr.GithubRef = 'https://github.com/chenchals/thomas-jpsth.git';

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
    [rho,pval] = corr(xMat,yMat,'type',corrMethod);
    [rho_pval] = [diag(rho),diag(pval)];
    n = size(xMat,1);
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
    if isstruct(dat)
      dat = struct2table(dat,'AsArray',true);
    end
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

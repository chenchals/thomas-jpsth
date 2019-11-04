function [spkCorr] = getSpkCorrForPair(cellPair,xSpkTimes,ySpkTimes,evntTimes,trialTypes,...
    conditions,alignNames,alignEvents,alignTimeWin,firstSortEvent)
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
    fx_mvsum = @(rasters,win) cellfun(@(x) movsum(double(x),win,2,'Endpoints','fill'),rasters,'UniformOutput',false);

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
            tempSpkCorr = table();
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
                alignTime = alignTime(selTrialsSorted);
                XAligned = SpikeUtils.alignSpikeTimes(units.(XCellId)(selTrialsSorted),alignTime, alignedTimeWin);
                YAligned = SpikeUtils.alignSpikeTimes(units.(YCellId)(selTrialsSorted),alignTime, alignedTimeWin);
                tempRast = SpikeUtils.rasters(XAligned,alignedTimeWin);
                XRasters = tempRast.rasters;
                tempRast = SpikeUtils.rasters(YAligned,alignedTimeWin);
                YRasters = tempRast.rasters;
                rasterBins = tempRast.rasterBins;
                % Z-Score using baseline mean and std
                
                clear tempRast;
                % gather vars fo corrSpk_PAIR_xxxx.mat
                opts(evId,1).condition = {condition};
                opts(evId,1).alignedName = {alignedName};
                opts(evId,1).alignedEvent = {alignedEvent};
                opts(evId,1).alignedTimeWin = {alignedTimeWin};
                opts(evId,1).alignTime = {alignTime}; 
                opts(evId,1).trialNosByCondition = selTrialsSorted;
                opts(evId,1).xCellSpikeTimes = {XAligned}; %#ok<*AGROW>
                opts(evId,1).yCellSpikeTimes = {YAligned};
                opts(evId,1).rasterBins = {rasterBins};
                opts(evId,1).xRasters = {XRasters};
                opts(evId,1).yRasters = {YRasters};
                opts(evId,1).firstSortByName = {firstSortName};
                opts(evId,1).firstSortByTime = selTrlsTblSorted.firstSortTime;
                %% Z scores
                % get Z score for given sample using means and stds specified for each trial
                fx_getZscore = @(rast,meanStd) arrayfun(@(x) (double(rast(x,:))- meanStd{1}(1))./meanStd{1}(2), (1:size(rast,1))','UniformOutput',false);
                opts(evId,1).xBaselineMeanStd = {units.X_baselineMeanStd};
                opts(evId,1).yBaselineMeanStd = {units.Y_baselineMeanStd};
                opts(evId,1).xTrialMeanStd = {units.X_trialMeanStd};
                opts(evId,1).yTrialMeanStd = {units.Y_trialMeanStd};
                opts(evId,1).xRasters_Z_baseline = fx_getZscore(XRasters,units.X_baselineMeanStd);
                opts(evId,1).xRasters_Z_trial = fx_getZscore(XRasters,units.X_trialMeanStd);
                opts(evId,1).yRasters_Z_baseline = fx_getZscore(YRasters,units.Y_baselineMeanStd);
                opts(evId,1).yRasters_Z_trial = fx_getZscore(YRasters,units.Y_trialMeanStd);
                
                %% compute spike corrs.....
                
                
            end % for alignEvents
            %tempSpkCorr.Properties.RowNames = alignNames;
            spkCorr.(condition) = [tempSpkCorr struct2table(opts,'AsArray',true)];
        catch mE
            getReport(mE)
            continue
        end
    end % for conditions
    %% Save for the current pair
    spkCorr.cellPairInfo = cellPair;
    
    
            out = struct();
        colNames = {'condition','alignedName','alignedEvent','alignedTimeWin',...
            'trialNosByCondition','alignTime','xCellSpikeTimes','yCellSpikeTimes',...
            'rasterBins','xRasters','yRasters','firstSortByName','firstSortByTime'...
            };

    
    
    
    
    
    
    
    spkCorr.GithubRef = 'https://github.com/chenchals/thomas-jpsth.git';

end


function [satJpsth] = getSatJpsth(jpsthPair,xSpkTimes,ySpkTimes,evntTimes,trialTypes,...
    conditions,alignNames,alignEvents,alignTimeWin,psthBinWidthMs,coincidenceBinWidthMs)

% jpsthPair : the pair row in the table JPSTH_PAIRS_CellInfoDB.mat
% xSpkTimes : X-axis cell: spike times as cell array of nTrials
% ySpkTimes : Y-axis cell: spike times as cell array of nTrials
% evntTimes : event times for session from TrialEventTimesDB.mat
% trialTypes : trial types for session from TrialTypesDB.mat
% conditions : Trial type names: AccurateCorrect, FastCorrect... from
%              TrialTypesDB
% alignNames : Name of the aligned event window : Baseline, Visual,
%              postSaccade
% alignEvents : Name of the event to align on ... from
%               TrialEventTimesDB.mat 
% alignTimeWin : A cell array of time windows for epochs corres[ponding to
%                align events
% psthBinWidthMs : binwidth for PSTH, Same bins will be used for JPSTH mat
% coincidenceBinWidthMs: binwidth for coincidence histogram.

warning('off');
% ignore processing if the sel. trials are below thisNum.
nTrialsThreshold = 10;

units = struct();
satJpsth = struct();

XCellId = ['DSP' jpsthPair.X_unit{1}];
YCellId = ['DSP' jpsthPair.Y_unit{1}];

units.(XCellId) = xSpkTimes;
units.(YCellId) = ySpkTimes;

%% For each condition
for cond = 1:numel(conditions)
    try
        % incase something breaks continue...
        condition = conditions{cond};
        selTrials = trialTypes.(condition){:};
        if isempty(selTrials)
            satJpsth.(condition) = [];
            continue;
        end
        
        %% Mutually Exclusive trila for Choice/Timing Errors
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
        if isempty(selTrials) || numel(selTrials) < nTrialsThreshold
            satJpsth.(condition) = [];
            continue;
        end
        
        %% check if trials need to be dropped due to poor_isolation..
        trRem = jpsthPair.X_trRemSAT{1};
        if ~isempty(trRem)
            selTrials(trRem(1):trRem(2)) = 0;
        end
        trRem = jpsthPair.Y_trRemSAT{1};
        if ~isempty(trRem)
            selTrials(trRem(1):trRem(2)) = 0;
        end
        if isempty(selTrials)
            satJpsth.(condition) = [];
            continue;
        end
        
        %% for each aligned event
        tempJpsth = table();
        opts = struct();
        for evId = 1:numel(alignEvents)
            alignedEvent = alignEvents{evId};
            alignedTimeWin = alignTimeWin{evId};
            alignedName = alignNames{evId};
            alignTime = evntTimes.CueOn{1};
            if ~strcmp(alignedEvent,'CueOn')
                alignTime = alignTime + evntTimes.(alignedEvent){1}(:);
            end
            alignTime = alignTime(selTrials);
            XAligned = SpikeUtils.alignSpikeTimes(units.(XCellId)(selTrials),alignTime, alignedTimeWin);
            YAligned = SpikeUtils.alignSpikeTimes(units.(YCellId)(selTrials),alignTime, alignedTimeWin);
            temp = SpikeUtils.jpsth(XAligned, YAligned, alignedTimeWin, psthBinWidthMs, coincidenceBinWidthMs);
            tempJpsth(evId,:) = struct2table(temp,'AsArray',true);
            %jer = SpikeUtils.jeromiahJpsth(XAligned, YAligned, alignedTimeWin, binWidth, coincidenceBins);
            opts(evId,1).xCellSpikeTimes = {XAligned}; %#ok<*AGROW>
            opts(evId,1).yCellSpikeTimes = {YAligned};
            opts(evId,1).trialNosByCondition = {selTrials};
            opts(evId,1).condition = {condition};
            opts(evId,1).alignedName = {alignedName};
            opts(evId,1).alignedEvent = {alignedEvent};
            opts(evId,1).alignedTimeWin = {alignedTimeWin};
            opts(evId,1).alignTime = {alignTime};
            opts(evId,1).binWidth = psthBinWidthMs;
            opts(evId,1).coincidenceBins = coincidenceBinWidthMs;
        end % for alignEvents
        tempJpsth.Properties.RowNames = alignNames;
        satJpsth.(condition) = [tempJpsth struct2table(opts,'AsArray',true)];
    catch mE
        disp(mE)
        continue
    end
end % for conditions
%% Save for the current pair
satJpsth.cellPairInfo = jpsthPair;
satJpsth.GithubRef = 'https://github.com/chenchals/thomas-jpsth.git';

end

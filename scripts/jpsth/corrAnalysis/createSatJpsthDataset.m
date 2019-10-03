function [] = createSatJpsthDataset(area1,area2)
% CREATESATJPSTHDATASET Create complete dataset for all the pairs matching the criteria for
% area1, area2
% Expects the following files in specified location:
% Inofrmation about all possible cell pairs:
%        'dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat'
% Trial types of all sessions (Accurate, Fast, Correct,...):
%        'dataProcessed/dataset/TrialTypesDB.mat'
% Event times for all trials and all sessions:
%         'dataProcessed/dataset/TrialEventTimesDB.mat'
% Get resptime and set it as SaccadePrimaryTempo
%         'dataProcessed/dataset/binfo_moves_SAT.mat'
% Spike time data for all units of all sessions:
%         'dataProcessed/dataset/spikes_SAT.mat'
% To generate pairs use:
%   for SEF-SEF pairs --> createSatJpsthDataset('SEF','FEF')
%   for SEF-SC pairs --> createSatJpsthDataset('SEF','SC')

%%
warning('off');
% ignore processing if the sel. trials are below thisNum.
nTrialsThreshold = 10;

%% Options for JPSTH computation
binWidth = 5;% use 1 ms for JPSTH computation
% -25 to +25 ms
coincidenceBins = 25/binWidth;

rootAnalysisDir = ['dataProcessed/analysis/JPSTH-' num2str(binWidth,'%dms')];
datasetDir = 'dataProcessed/dataset';
jpsthResultsDir = fullfile(rootAnalysisDir,['jpsth_' area1 '-' area2]);
if ~exist(jpsthResultsDir, 'dir')
    mkdir(jpsthResultsDir);
end

% Info files
jpshPairsFile = fullfile(datasetDir,'JPSTH_PAIRS_CellInfoDB.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');
% to get resptime and set it as SaccadePrimaryTempo
binfoFile = fullfile(datasetDir,'binfo_moves_SAT.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignNames = {'Baseline','Visual', 'PostSaccade'};
alignEvents = {'CueOn', 'CueOn','SaccadePrimaryTempo'};
alignTimeWin = {[-700 100], [-100 400], [-100 400]};

%conditions
% (Accurate|Fast)*(Correct|ErrorHold|ErrorChoice|ErrorTiming|ErrorNoSaccade)
availConditions = {'Accurate'; 'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
    'Fast';     'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'};
%% Load Variables
% load variable: JpsthPairsCellInfo
jpsthCellPairs = load(jpshPairsFile);
jpsthCellPairs = jpsthCellPairs.JpsthPairCellInfoDB;
% load all spike times
spikeTimes = load(spikeTimesFile);
spikeTimes = spikeTimes.spikes;
% load TrialTypes
trialTypes = load(trialTypesFile);
trialTypes = trialTypes.TrialTypesDB;
% load TrialEventTimes
trialEventTimes = load(trialEventTimesFile);
trialEventTimes = trialEventTimes.TrialEventTimesDB;
% The variable SaccadePrimaryTempo is NaN, so replace with one from binfo
binfo = load(binfoFile);
temp = table();
temp.session = {binfo.binfo.SAT.session}';
temp.resptime = {binfo.binfo.SAT.resptime}';
temp = innerjoin(trialEventTimes,temp);
trialEventTimes.SaccadePrimaryTempo = temp.resptime;
clearvars temp binfo

%% Filter cell pairs for the areas of interest
jpsthCellPairs = jpsthCellPairs(((contains(jpsthCellPairs.X_area,area1) & contains(jpsthCellPairs.Y_area,area2)) ...
    | (contains(jpsthCellPairs.X_area,area2) & contains(jpsthCellPairs.Y_area,area1))),:);
sessions = unique(jpsthCellPairs.X_sess);
assert(isequal(sessions,unique(jpsthCellPairs.Y_sess)),'********Fatal: Error X-xell sessions and Y-cell sessions do not match');

%% Filter Trial Types and Trial Event Times for the sessions of interest above
trialTypes = trialTypes(contains(trialTypes.session,sessions),:);
trialEventTimes = trialEventTimes(contains(trialEventTimes.session,sessions),:);

%% For each JPSH cell pair do JPSTH
% see doc pctRunOnAll
pctRunOnAll warning off;
parfor s = 1:size(jpsthCellPairs,1)
    units = struct();
    jpsthPair = jpsthCellPairs(s,:);
    sessionTrialEventTimes = trialEventTimes(contains(trialEventTimes.session,jpsthPair.X_sess),:);
    sessionTrialTypes = trialTypes(contains(trialTypes.session,jpsthPair.X_sess),:);
    
    tempConditions = struct();
    pairUid = jpsthPair.Pair_UID{1};
    XCellId = ['DSP' jpsthPair.X_unit{1}];
    YCellId = ['DSP' jpsthPair.Y_unit{1}];
    pairFilename = ['JPSTH-' pairUid];
    units.(XCellId) = spikeTimes(jpsthPair.X_unitNum).SAT';
    units.(YCellId) = spikeTimes(jpsthPair.Y_unitNum).SAT';
    
    %% For each condition
    for cond = 1:numel(availConditions)
        try
            % incase something breaks continue...
            condition = availConditions{cond};
            selTrials = sessionTrialTypes.(condition){:};
            if isempty(selTrials)
                tempConditions.(condition) = [];
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
                selTrials(sessionTrialTypes.(otherCondition){:}) = 0;
            end
            if isempty(selTrials) || numel(selTrials) < nTrialsThreshold
                tempConditions.(condition) = [];
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
                tempConditions.(condition) = [];
                continue;
            end
            
            %% for each aligned event
            tempJpsth = table();
            opts = struct();
            for evId = 1:numel(alignEvents)
                alignedEvent = alignEvents{evId};
                alignedTimeWin = alignTimeWin{evId};
                alignedName = alignNames{evId};
                alignTime = sessionTrialEventTimes.CueOn{1};
                if ~strcmp(alignedEvent,'CueOn')
                    alignTime = alignTime + sessionTrialEventTimes.(alignedEvent){1}(:);
                end
                alignTime = alignTime(selTrials);
                XAligned = SpikeUtils.alignSpikeTimes(units.(XCellId)(selTrials,:),alignTime, alignedTimeWin);
                YAligned = SpikeUtils.alignSpikeTimes(units.(YCellId)(selTrials,:),alignTime, alignedTimeWin);
                temp = SpikeUtils.jpsth(XAligned, YAligned, alignedTimeWin, binWidth, coincidenceBins);
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
                opts(evId,1).binWidth = binWidth;
                opts(evId,1).coincidenceBins = coincidenceBins;
            end % for alignEvents
            tempJpsth.Properties.RowNames = alignNames;
            tempConditions.(condition) = [tempJpsth struct2table(opts,'AsArray',true)];
        catch mE
            disp(mE)
            continue
        end
    end % for conditions
    %% Save for the current pair
    tempConditions.cellPairInfo = jpsthPair;
    tempConditions.GithubRef = 'https://github.com/chenchals/thomas-jpsth.git';
    oFn = fullfile(jpsthResultsDir,[pairFilename '.mat']);
    fprintf('Saving processed pair : %s\n',oFn);
    saveJpsthData(oFn,tempConditions);
    %save(oFn,'-v7.3','-struct','tempConditions');
    
end

end

function [] = saveJpsthData(oFn,varToSave)
tempConditions = varToSave;
save(oFn,'-v7.3','-struct','tempConditions');
end

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

%% Options for JPSTH computation
psthBinWidthMs = 10;
% -25 to +25 ms
coincidenceBinWidthMs = 25;

rootAnalysisDir = ['dataProcessedL/analysis/JPSTH-' num2str(psthBinWidthMs,'%dms')];
datasetDir = 'dataProcessed/dataset';
jpsthResultsDir = fullfile(rootAnalysisDir,['jpsth_' area1 '-' area2],'mat');
if ~exist(jpsthResultsDir, 'dir')
    mkdir(jpsthResultsDir);
end

%% Files for getting data to run JPSTH
jpshPairsFile = fullfile(datasetDir,'JPSTH_PAIRS_CellInfoDB.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');

%% Load data variable: JpsthPairsCellInfo
jpsthCellPairs = load(jpshPairsFile);
jpsthCellPairs = jpsthCellPairs.JpsthPairCellInfoDB;
% Load data variable: spike times
spikesSat = load(spikeTimesFile);
spikesSat = {spikesSat.spikes.SAT};
% Load data variable: TrialTypes 
sessionTrialTypes = load(trialTypesFile);
sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
% Load data variable: TrialEventTimes 
sessionEventTimes = load(trialEventTimesFile);
sessionEventTimes = sessionEventTimes.TrialEventTimesDB;

%% Filter cell pairs for the areas of interest
jpsthCellPairs = jpsthCellPairs(...
    ((contains(jpsthCellPairs.X_area,area1) & contains(jpsthCellPairs.Y_area,area2)) ...
   | (contains(jpsthCellPairs.X_area,area2) & contains(jpsthCellPairs.Y_area,area1))),...
     :);
sessions = unique(jpsthCellPairs.X_sess);
assert(isequal(sessions,unique(jpsthCellPairs.Y_sess)),'********Fatal: Error X-xell sessions and Y-cell sessions do not match');

%% alignment:
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignNames = {'Baseline','Visual','PostSaccade','PostReward'};
alignEvents = {'CueOn','CueOn','SaccadePrimary','RewardTime'};
alignTimeWin = {[-600 100],[-200 400],[-100 500],[-200 700]};
firstSortEvent = {'SaccadePrimary','SaccadePrimary','SaccadeSecond',[]};
%%conditions
% (Accurate|Fast)*(Correct|ErrorHold|ErrorChoice|ErrorTiming|ErrorNoSaccade)
conditionsTbl = table();
conditionsTbl.conditions = {
     'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
     'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
     };
% sort trials first event
conditionsTbl.trialsSortFirst = {
    'SaccadePrimary';'SaccadePrimary';'SaccadePrimary';
    'SaccadePrimary';'SaccadePrimary';'SaccadePrimary';
    };
% sort trials next event
conditionsTbl.trialsSortNext = {
    'RewardTime';'SaccadeSecond';[];
    'RewardTime';'SaccadeSecond';[];
    };

%% For each JPSH cell pair do JPSTH
% see doc pctRunOnAll
pctRunOnAll warning off;
nPairs = size(jpsthCellPairs,1);
parfor p = 1:nPairs
    jpsthPair = jpsthCellPairs(p,:); %#ok<*PFBNS>
    sess = jpsthPair.X_sess{1};
    % must be cell array of ntrials by 1
    xSpkTimes = spikesSat{jpsthPair.X_unitNum}';
    ySpkTimes = spikesSat{jpsthPair.Y_unitNum}';
    evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
    trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
    satJpsth =  getSatJpsthForPair(jpsthPair,xSpkTimes,ySpkTimes,...
                            evntTimes,trialTypes,conditionsTbl.conditions,...
                            alignNames,alignEvents,alignTimeWin,firstSortEvent,...
                            psthBinWidthMs,coincidenceBinWidthMs);
    oFn = fullfile(jpsthResultsDir,['JPSTH-' jpsthPair.Pair_UID{1} '.mat']);
    saveJpsthData(oFn,satJpsth);
end

end

function [] = saveJpsthData(oFn,varToSave)
fprintf('Saving file : %s\n',oFn)
tempConditions = varToSave;
save(oFn,'-v7.3','-struct','tempConditions');
end


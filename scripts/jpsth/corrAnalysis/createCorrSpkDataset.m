function [] = createCorrSpkDataset(area1,area2,wavDir)
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
monkIdsToDo = {'D','E'};

%% Options for Spk Corr computation
rootAnalysisDir = 'dataProcessed/analysis/spkCorr';
datasetDir = 'dataProcessed/dataset';
resultsDir = fullfile(rootAnalysisDir,['spkCorr_' area1 '-' area2],'mat');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Files for getting data to run JPSTH
pairsFile = fullfile(datasetDir,'JPSTH_PAIRS_CellInfoDB.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');

%% Load data variable: JpsthPairsCellInfo
cellPairs = load(pairsFile);
cellPairs = cellPairs.JpsthPairCellInfoDB;
cellPairs = cellPairs(ismember([cellPairs.X_monkey],monkIdsToDo),:);
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
cellPairs = cellPairs(...
    ((strcmp(cellPairs.X_area,area1) & strcmp(cellPairs.Y_area,area2)) ...
   | (strcmp(cellPairs.X_area,area2) & strcmp(cellPairs.Y_area,area1))),...
     :);
sessions = unique(cellPairs.X_sess);
assert(isequal(sessions,unique(cellPairs.Y_sess)),'********Fatal: Error X-Unit sessions and Y-Unit sessions do not match');

%% alignment:
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignNames = {'Baseline','Visual','PostSaccade','PostReward'};
alignEvents = {'CueOn','CueOn','SaccadePrimary','RewardTime'};
alignTimeWin = {[-600 100],[-200 400],[-100 500],[-200 700]};
% first sort by accurate/fast?
firstSortEvent = {'SaccadePrimary','SaccadePrimary','SaccadeSecond',[]};
%%conditions
conditionsTbl = table();
conditionsTbl.conditions = {
     'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
     'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
     };
% % sort trials first event
% conditionsTbl.trialsSortFirst = {
%     'SaccadePrimary';'SaccadePrimary';'SaccadePrimary';
%     'SaccadePrimary';'SaccadePrimary';'SaccadePrimary';
%     };
% % sort trials next event
% conditionsTbl.trialsSortNext = {
%     'RewardTime';'SaccadeSecond';[];
%     'RewardTime';'SaccadeSecond';[];
%     };

%% For each JPSH cell pair do JPSTH
% see doc pctRunOnAll
pctRunOnAll warning off;
nPairs = size(cellPairs,1);
for p = 1:nPairs
    cellPair = cellPairs(p,:); %#ok<*PFBNS>
    sess = cellPair.X_sess{1};
    % must be cell array of ntrials by 1
    xSpkTimes = spikesSat{cellPair.X_unitNum}';
    ySpkTimes = spikesSat{cellPair.Y_unitNum}';
    evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
    trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
    spkCorr = getSpkCorrForPair(cellPair,xSpkTimes,ySpkTimes,...
                            evntTimes,trialTypes,conditionsTbl.conditions,...
                            alignNames,alignEvents,alignTimeWin,firstSortEvent,...
                            wavDir);    
    oFn = fullfile(resultsDir,['spkCorr_' cellPair.Pair_UID{1} '.mat']);
    saveSpkCorrData(oFn,spkCorr);
end

end

function [] = saveSpkCorrData(oFn,varToSave)
    fprintf('Saving file : %s\n',oFn)
    tempConditions = varToSave;
    save(oFn,'-v7.3','-struct','tempConditions');
end


% Get sat SDFs for a given unit

spkCorrDir = 'dataProcessed/analysis/11-18-2019/spkCorr';
% for loading event times and trial types
datasetDir = 'dataProcessed/dataset';
% Files from which to load spik corrs and compute SDFs
spkCorrFile = fullfile(spkCorrDir,'summary/spkCorrAllPairsStaticNew.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
unitInfoStatsFile = fullfile(datasetDir,'dataNeurophys_SAT.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');
%% load once in a script
if ~exist('spkCorr','var')
    % Load data variable: spike correlations
    spkCorr = load(spkCorrFile);
    % Load data variable: spike times
    spikesSat = load(spikeTimesFile);
    spikesSat = spikesSat.spikesSAT;
    % Load unit Info for knowing which trials to be removed of any
    unitInfoAll = load(unitInfoStatsFile,'unitInfo');
    unitInfoAll = unitInfoAll.unitInfo;
    % Load data variable: TrialTypes
    sessionTrialTypes = load(trialTypesFile);
    sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
    % Load data variable: TrialEventTimes
    sessionEventTimes = load(trialEventTimesFile);
    sessionEventTimes = sessionEventTimes.TrialEventTimesDB;
end
%% setup vars
     useAlignment.events = {'CueOn','SaccadePrimary','RewardTime'};
     useAlignment.timeWins = {[-600 400],[-200 600],[-100 400]};
     useAlignment.names = {'Visual','PostSaccade','PostReward'};
     useAlignment.sortEventNames = {'SaccadePrimary','SaccadePrimary','SaccadePrimary'};

     unitNums = [13,14,26,10];
     
     
     for ii = 1:numel(unitNums)
         unitNum = unitNums(ii);
         unitInfo = unitInfoAll(unitNum,:);
         sess = unitInfo.sess{1};
         
         useUnit.unitNum = unitNum;
         useUnit.spkTimes = spikesSat{unitNum}';
         evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
         useTrials.trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
         useTrials.trialsToRemove = unitInfo.trRemSAT{1};
         useTrials.minTrialCount = 1;
         [outTbl] = getUnitSatSdf(useUnit,evntTimes,useTrials,useAlignment);
         
     end

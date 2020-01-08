function [] = createSatSdfsAll()
%% Create a dataset of unit SDFs for *all* units by SAT condition
% use the call to [outTbl] = getUnitSatSdf(useUnit,evntTimes,useTrials,useAlignment)
% Requires the following files for getting events, spikes, etc
% Use data from:
% dataset/spikes_SAT.mat
% dataset/dataNeurophys_SAT.mat
% dataset/TrialTypesDB.mat
% dataset/TrialEventTimesDB.mat
% see also GETUNITSATSDF
%
    %% Load data needed for computing SDFs
    datasetDir = 'dataProcessed/dataset';
    satSdfDir = 'dataProcessed/dataset/satSdfs';
    if ~exist(satSdfDir,'dir')
        mkdir(satSdfDir);
    end
    % Load data variable: spike times
    spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
    spikesSat = load(spikeTimesFile);
    spikesSat = spikesSat.spikesSAT;
    % Load unit Info for knowing which trials to be removed of any
    unitInfoStatsFile = fullfile(datasetDir,'dataNeurophys_SAT.mat');
    unitInfoAll = load(unitInfoStatsFile,'unitInfo');
    unitInfoAll = unitInfoAll.unitInfo;
    % Load data variable: TrialTypes
    trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
    sessionTrialTypes = load(trialTypesFile);
    sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
    % Load data variable: TrialEventTimes
    trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');
    sessionEventTimes = load(trialEventTimesFile);
    sessionEventTimes = sessionEventTimes.TrialEventTimesDB;
    fprintf('loaded data variables...\n');
    %% Compute SDFs by SAT condition for each unit
    unitNums = unitInfoAll.unitNum;
    % argument: alignmentStruct
    useAlignment.events = {'CueOn','SaccadePrimary','RewardTime'};
    useAlignment.timeWins = {[-600 400],[-200 600],[-100 400]};
    useAlignment.names = {'Visual','PostSaccade','PostReward'};
    useAlignment.sortEventNames = {'SaccadePrimary','SaccadePrimary','SaccadePrimary'};
    %%
    parfor uu = 1:numel(unitNums)
        useTrials = struct();
        useUnit = struct();
        unitNum = unitNums(uu);
        sess = unitInfoAll.sess{unitInfoAll.unitNum == unitNum}; %#ok<*PFBNS>
        % only if this session has eventTimes and trialTypes process unit
        if any(ismember(sessionEventTimes.session,sess)) ...
                && any(ismember(sessionTrialTypes.session,sess))
            % argument eventTimes
            useEventTimes = sessionEventTimes(ismember(sessionEventTimes.session,sess),:);
            % argument: trialsStruct
            useTrials.trialTypes = sessionTrialTypes(ismember(sessionTrialTypes.session,sess),:);
            useTrials.trialsToRemove = unitInfoAll.trRemSAT{unitInfoAll.unitNum == unitNum};
            useTrials.minTrialCount = 1;
            % argument: unitStruct
            useUnit.unitNum = unitNum;
            useUnit.spkTimes = spikesSat{unitNum}';
            sdfTbl = getUnitSatSdf(useUnit,useEventTimes,useTrials,useAlignment);
            oFn = fullfile(satSdfDir,sprintf('Unit_%03d.mat',unitNum));
            saveSatSdf(oFn,sdfTbl,useAlignment,useTrials,useEventTimes);
        end
    end
end

function [] = saveSatSdf(oFn,sdfs,alignment,trialTypes,eventTimes)
  temp.sdfs = sdfs;
  temp.alignment = alignment;
  temp.trialTypes = trialTypes;
  temp.eventTimes = eventTimes;
  save(oFn,'-v7.3','-struct','temp');
end

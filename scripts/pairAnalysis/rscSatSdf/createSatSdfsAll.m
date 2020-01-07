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
    mFilename = mfilename;
    datasetDir = 'dataProcessed/dataset';
    satSdfFilename = 'dataProcessed/dataset/SDFs_SAT.mat';
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
    %% Compute SDFs by SAT condition for each unit
    unitNums = unitInfoAll.unitNum;
    % argument: alignmentStruct
    useAlignment.events = {'CueOn','SaccadePrimary','RewardTime'};
    useAlignment.timeWins = {[-600 400],[-200 600],[-100 400]};
    useAlignment.names = {'Visual','PostSaccade','PostReward'};
    useAlignment.sortEventNames = {'SaccadePrimary','SaccadePrimary','SaccadePrimary'};
    %%
    if ~exist(satSdfFilename,'file')
        save(satSdfFilename,'-v7.3','useAlignment','mFilename')
    end
    for uu = 1:numel(unitNums)
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
            [outTbl] = getUnitSatSdf(useUnit,useEventTimes,useTrials,useAlignment);
            appendVar(satSdfFilename,outTbl,unitNum);
        end
    end 
end

function [] = appendVar(oFn,sdfTbl,unitNum)
   temp.(sprintf('Unit_%03d',unitNum)) = sdfTbl;
   save(oFn,'-v7.3','-append','-struct','temp')
end


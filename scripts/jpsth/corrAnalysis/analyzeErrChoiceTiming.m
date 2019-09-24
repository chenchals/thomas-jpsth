
JpsthPairCellInfoDB = load('dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat');
JpsthPairCellInfoDB=JpsthPairCellInfoDB.JpsthPairCellInfoDB;
% SEF error units
errUnits = load('dataProcessed/dataset/ID_Signal_Error.mat');
choiceErrUnits = errUnits.unit_Signal_ChoiceError;
timingErrUnits = errUnits.unit_Signal_TimingError;

choiceErrPairsIdx = arrayfun(@(x) find((JpsthPairCellInfoDB.X_unitNum==x|JpsthPairCellInfoDB.Y_unitNum==x)...
    & strcmp(JpsthPairCellInfoDB.X_area,'SEF') ...
    & strcmp(JpsthPairCellInfoDB.Y_area,'SEF')),[choiceErrUnits.unitNum],'UniformOutput',false);
choiceErrPairsIdx = cell2mat(choiceErrPairsIdx);
timingErrPairsIdx = arrayfun(@(x) find((JpsthPairCellInfoDB.X_unitNum==x|JpsthPairCellInfoDB.Y_unitNum==x)...
    & strcmp(JpsthPairCellInfoDB.X_area,'SEF') ...
    & strcmp(JpsthPairCellInfoDB.Y_area,'SEF')),[timingErrUnits.unitNum],'UniformOutput',false);
timingErrPairsIdx = cell2mat(timingErrPairsIdx);

% to create JPSTH data we need
% spikes
allSpikes = load('dataProcessed/dataset/spikes_SAT.mat');
spikesSat = {allSpikes.spikes.SAT};
% event times
sessionEventTimes = load('dataProcessed/dataset/TrialEventTimesDB.mat');
sessionEventTimes = sessionEventTimes.TrialEventTimesDB;
% trial types: accurate, fast, accurateErr
sessionTrialTypes = load('dataProcessed/dataset/TrialTypesDB.mat');
sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
% conditions:
conditions = {'Accurate'; 'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
    'Fast';     'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'};
% alignment:
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignNames = {'Baseline','Visual', 'PostSaccade'};
alignEvents = {'CueOn', 'CueOn','SaccadePrimaryTempo'};
alignTimeWin = {[-700 100], [-100 500], [-100 600]};

psthBinWidthMs = 5;
coincidenceBinWidthMs = 25;

if ~exist('temp','dir')
    mkdir temp
end
%% process for pairs
errPairs = {choiceErrPairsIdx, timingErrPairsIdx};
pairDirs = {'dataProcessed/analysis/SEF-PAPER/CHOICE_ERR_PAIRS'
    'dataProcessed/analysis/SEF-PAPER/TIMING_ERR_PAIRS'
    };
for jj = 1:numel(errPairs)
    errPair = errPairs{jj};
    pairDir = pairDirs{jj};
    if ~exist(pairDir,'dir')
        mkdir(pairDir)
    end
    parfor p = 1:1 % numel(choiceErrPairsIdx)
        jpsthPair = JpsthPairCellInfoDB(choiceErrPairsIdx(p),:); %#ok<*PFBNS>
        sess = jpsthPair.X_sess{1};
        % must be cell array of ntrials by 1
        xSpkTimes = spikesSat{jpsthPair.X_unitNum}';
        ySpkTimes = spikesSat{jpsthPair.Y_unitNum}';
        evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
        trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
        satJpsth =  getSatJpsth(jpsthPair,xSpkTimes,ySpkTimes,evntTimes,trialTypes,conditions,alignNames,alignEvents,alignTimeWin,psthBinWidthMs,coincidenceBinWidthMs);
        oFn = fullfile(pairDir,['mat/JPSTH-' jpsthPair.Pair_UID '.mat']);
        saveJpsthData(oFn,satJpsth);
    end
    
end

function [] = saveJpsthData(oFn,varToSave)
    fprintf('Saving processed pair : %s\n',oFn);
    tempConditions = varToSave;
    save(oFn,'-v7.3','-struct','tempConditions');
end


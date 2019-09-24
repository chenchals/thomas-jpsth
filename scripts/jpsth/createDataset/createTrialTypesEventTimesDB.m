% Needs 'JPSTH_PAIRS_CellInfoDB.mat' created by createJpsthPairCellInfoDB()
% see also: CREATEJPSTHPAIRCELLINFODB

datasetDir = 'dataProcessed/dataset';
jpsthPairsDBFile = fullfile(datasetDir,'JPSTH_PAIRS_CellInfoDB.mat');
% Output dataset files
trialTypesFile = fullfile(datasetDir,'TrialTypesDB_2.mat');
TrialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB_2.mat');

%% Process for trial event times and trial types

jpsthPairsDB = load(jpsthPairsDBFile);
jpsthPairsDB = jpsthPairsDB.JpsthPairCellInfoDB;

sessionMatFiles = unique(jpsthPairsDB(:,{'X_sess','matDatafile'}));

% Create TrialTypes for each session
% Regexp for all vars ending in '_' and not starting with Eye or Pupil
vars2LoadRegEx = '.*_$(?<!^(Eye|Pupil).*)|saccLoc|SRT';
TrialTypesDB = struct();
TrialEventTimesDB = struct();
nSess = size(sessionMatFiles,1);
for s = 1:nSess
    sessionFile = sessionMatFiles.matDatafile{s};
    sessName = sessionMatFiles.X_sess{s};
    vars = load(sessionFile,'-regexp',vars2LoadRegEx);
    
    fprintf('Doing session [%i] of [%i]: [%s]...\n',s,nSess,sessionFile);    
 
    nTrials = size(vars.Correct_,1);
    TrialTypesDB.session{s,1} = sessName;
    %% Trial type and conditions
    % From Thomas' code
    % info(kk).condition = transpose(uint8(SAT_(:,1))); %1==accurate, 3==fast
    accurate = vars.SAT_(:,1) == 1;
    fast = vars.SAT_(:,1) == 3;
    correct = vars.Correct_(:,2) ==1;
    % From Thomas' code
    % Response information
    nosacc = vars.Errors_(:,2) == 1;
    err_hold = vars.Errors_(:,4) == 1;
    err_dir = vars.Errors_(:,5) == 1;
    err_time = vars.Errors_(:,6) == 1 | vars.Errors_(:,7) == 1;
    % Different Trial types
    TrialTypesDB.Accurate{s,1} = accurate;
    TrialTypesDB.AccurateCorrect{s,1} = accurate & correct;
    TrialTypesDB.AccurateErrorHold{s,1} = accurate & err_hold;
    TrialTypesDB.AccurateErrorChoice{s,1} = accurate & err_dir;
    TrialTypesDB.AccurateErrorTiming{s,1} = accurate & err_time;
    TrialTypesDB.AccurateErrorNoSaccade{s,1} = accurate & nosacc;
    TrialTypesDB.Fast{s,1} = fast;
    TrialTypesDB.FastCorrect{s,1} = fast & correct;
    TrialTypesDB.FastErrorHold{s,1} = fast & err_hold;
    TrialTypesDB.FastErrorChoice{s,1} = fast & err_dir;
    TrialTypesDB.FastErrorTiming{s,1} = fast & err_time;
    TrialTypesDB.FastErrorNoSaccade{s,1} = fast & nosacc;
    % Stimulus/Response Location
    TrialTypesDB.SingletonLoc{s,1} = vars.Target_(:,2);
    TrialTypesDB.ResponseLoc{s,1} = vars.saccLoc;
    
    %% SAT event times
    TrialEventTimesDB.session{s,1} = sessName;
    TrialEventTimesDB.TrialStart{s,1} = nan(nTrials,1);
    if isfield(vars,'TrialStart_')
        TrialEventTimesDB.TrialStart{s,1} = vars.TrialStart_(:,1);       
    end
    TrialEventTimesDB.CueOn{s,1} = nan(nTrials,1);
    if isfield(vars,'Target_')
        TrialEventTimesDB.CueOn{s,1} = vars.Target_(:,1);       
    end
    TrialEventTimesDB.FixAcquisition{s,1} = nan(nTrials,1);
    if isfield(vars,'FixAcqTime_')
        TrialEventTimesDB.FixAcquisition{s,1} = vars.FixAcqTime_(:,1);       
    end
    TrialEventTimesDB.TargetDeadline{s,1} = nan(nTrials,1);
    if isfield(vars,'SAT_')
        temp = vars.SAT_(:,3);
        temp(temp > 1000) = NaN;
        TrialEventTimesDB.TargetDeadline{s,1} = temp; 
        clearvars temp;
    end    
    TrialEventTimesDB.SaccadePrimaryTempo{s,1} = nan(nTrials,1);
    if isfield(vars,'SRT')
        TrialEventTimesDB.SaccadePrimaryTempo{s,1} = vars.SRT(:,1);       
    end     
    TrialEventTimesDB.ToneOn{s,1} = nan(nTrials,1);
    if isfield(vars,'ToneOn_')
        TrialEventTimesDB.ToneOn{s,1} = vars.ToneOn_(:,1);       
    end
    TrialEventTimesDB.RewardOn{s,1} = nan(nTrials,1);
    if isfield(vars,'JuiceOn_')
        TrialEventTimesDB.RewardOn{s,1} = vars.JuiceOn_(:,1);       
    end
end
TrialTypesDB = struct2table(TrialTypesDB);
TrialEventTimesDB = struct2table(TrialEventTimesDB);

save(trialTypesFile,'TrialTypesDB');
save(TrialEventTimesFile, 'TrialEventTimesDB');


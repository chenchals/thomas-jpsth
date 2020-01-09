% Needs 'JPSTH_PAIRS_CellInfoDB.mat' created by createJpsthPairCellInfoDB()
%       'dataBehavior_SAT.mat' created by Thomas (10/29/2019)
% see also: CREATEJPSTHPAIRCELLINFODB
% Modifications:
% 10/29/2019 : Using new file for dataBehabior_SAT.mat instead of
%              binfo_moves_SAT.mat 

datasetDir = 'dataProcessed/dataset';
jpsthPairsDBFile = fullfile(datasetDir,'JPSTH_PAIRS_CellInfoDB.mat');
binfoMovesFile = fullfile(datasetDir,'dataBehavior_SAT.mat');
% Output dataset files
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
TrialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');

%% Process for trial event times and trial types

jpsthPairsDB = load(jpsthPairsDBFile);
jpsthPairsDB = jpsthPairsDB.JpsthPairCellInfoDB;
sessionMatFiles = unique(jpsthPairsDB(:,{'X_sess','matDatafile'}));

% load binfo_moves_SAT.mat for
% 1. saccadePrimary, saccadeSecondary, and rewardOn 
temp = load(binfoMovesFile);
sessionNames = temp.binfoSAT.session;
sessionRewardT = temp.binfoSAT.rewtime;
sessionPrimarySaccadeT = temp.primarySaccade.resptime;
sessionSecondSaccadeT = temp.secondSaccade.resptime;

% Create TrialTypes for each session
% Regexp for all vars ending in '_' and not starting with Eye or Pupil
vars2LoadRegEx = '.*_$(?<!^(Eye|Pupil).*)|saccLoc|SRT';

nSess = size(sessionMatFiles,1);
out = struct();

parfor s = 1:nSess
    sessionFile = sessionMatFiles.matDatafile{s};
    sessName = sessionMatFiles.X_sess{s};
    vars = load(sessionFile,'-regexp',vars2LoadRegEx);
    fprintf('Doing session [%i] of [%i]: [%s]...\n',s,nSess,sessionFile);    
    
    nTrials = size(vars.Correct_,1);
    out(s).TrialTypesDB.session = sessName;
    %% Trial type and conditions
    % From Thomas' code
    % info(kk).condition = transpose(uint8(SAT_(:,1))); %1==accurate,
    % 3==fast, 4 = NaN
    nanTrial = vars.SAT_(:,1) == 4;
    accurate = vars.SAT_(:,1) == 1;
    fast = vars.SAT_(:,1) == 3;
    correct = vars.Correct_(:,2) ==1;
    % Do block number for each trial...
    blkStartIdx = find([0;diff(accurate-fast)]~=0);
    blkNumbers = (1:numel(blkStartIdx))';
    % session may not start with accurate or fast SAT_(:,1) ~= [1|3]
    if blkStartIdx(1)~=1
        blkStartIdx =[1;blkStartIdx]; %#ok<*AGROW>
        blkNumbers = [NaN;blkNumbers];
    end
    blkEndIdx = [blkStartIdx(2:end)-1;numel(accurate)];    
    trlBlkNumber = arrayfun(@(s,e,n) repmat(n,e-s+1,1),blkStartIdx,blkEndIdx,blkNumbers,'UniformOutput',false);
    trlBlkNumber = cell2mat(trlBlkNumber);   
     
    % From Thomas' code
    % Response information
    nosacc = vars.Errors_(:,2) == 1;
    err_hold = vars.Errors_(:,4) == 1;
    err_dir = vars.Errors_(:,5) == 1;
    err_time = vars.Errors_(:,6) == 1 | vars.Errors_(:,7) == 1;
    % Different Trial types
    
    out(s).TrialTypesDB.TrialBlockNum = trlBlkNumber;
    out(s).TrialTypesDB.Accurate = double(accurate);
    out(s).TrialTypesDB.AccurateCorrect = accurate & correct;
    out(s).TrialTypesDB.AccurateErrorHold= accurate & err_hold;
    out(s).TrialTypesDB.AccurateErrorChoice = accurate & err_dir;
    out(s).TrialTypesDB.AccurateErrorTiming = accurate & err_time;
    out(s).TrialTypesDB.AccurateErrorNoSaccade = accurate & nosacc;
    out(s).TrialTypesDB.Fast = double(fast);
    out(s).TrialTypesDB.FastCorrect = fast & correct;
    out(s).TrialTypesDB.FastErrorHold = fast & err_hold;
    out(s).TrialTypesDB.FastErrorChoice = fast & err_dir;
    out(s).TrialTypesDB.FastErrorTiming = fast & err_time;
    out(s).TrialTypesDB.FastErrorNoSaccade = fast & nosacc;
    %
    out(s).TrialTypesDB.Accurate(nanTrial) = NaN;
    out(s).TrialTypesDB.Fast(nanTrial) = NaN;
    
    % Stimulus/Response Location
    out(s).TrialTypesDB.SingletonLoc = vars.Target_(:,2);
    out(s).TrialTypesDB.ResponseLoc = vars.saccLoc;
    
    %% SAT event times
    out(s).TrialEventTimesDB.session = sessName;
    out(s).TrialEventTimesDB.TrialStart = nan(nTrials,1);
    out(s).TrialEventTimesDB.TrialBlockNum = trlBlkNumber;
    if isfield(vars,'TrialStart_')
        out(s).TrialEventTimesDB.TrialStart = vars.TrialStart_(:,1);       
    end
    out(s).TrialEventTimesDB.CueOn = nan(nTrials,1);
    if isfield(vars,'Target_')
        out(s).TrialEventTimesDB.CueOn = vars.Target_(:,1);       
    end
    out(s).TrialEventTimesDB.FixAcquisition = nan(nTrials,1);
    if isfield(vars,'FixAcqTime_')
        out(s).TrialEventTimesDB.FixAcquisition = vars.FixAcqTime_(:,1);       
    end
    out(s).TrialEventTimesDB.TargetDeadline = nan(nTrials,1);
    if isfield(vars,'SAT_')
        temp = vars.SAT_(:,3);
        temp(temp > 1000) = NaN;
        out(s).TrialEventTimesDB.TargetDeadline = temp; 
    end  
    out(s).TrialEventTimesDB.BellOn = nan(nTrials,1);
    if isfield(vars,'BellOn_')
        out(s).TrialEventTimesDB.BellOn = vars.BellOn_(:,1);       
    end
    out(s).TrialEventTimesDB.JuiceOn = nan(nTrials,1);
    if isfield(vars,'JuiceOn_')
        out(s).TrialEventTimesDB.JuiceOn = vars.JuiceOn_(:,1);       
    end
    % get from binfo and movespp :
    %   primary saccde, second saccade, and reward time 
    idx = strcmp(sessionNames,sessName);
    out(s).TrialEventTimesDB.SaccadePrimary = sessionPrimarySaccadeT{idx}';
    out(s).TrialEventTimesDB.SaccadeSecond = sessionSecondSaccadeT{idx}';
    out(s).TrialEventTimesDB.RewardTime = (sessionPrimarySaccadeT{idx} + sessionRewardT{idx})';
    
end
TrialTypesDB = struct2table([out.TrialTypesDB],'AsArray',true);
TrialEventTimesDB = struct2table([out.TrialEventTimesDB],'AsArray',true);

save(trialTypesFile,'TrialTypesDB');
save(TrialEventTimesFile, 'TrialEventTimesDB');


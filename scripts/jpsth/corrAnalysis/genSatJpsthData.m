function [] = genSatJpsthData(area1,area2)
% use tic;genJpsthData('SEF','FEF');toc;tic;genJpsthData('SEF','SC');toc 
% to generate all pairs
%% Cross-area JPSTH analysis of cell pairs
% SEF_FEF 
% [includes V,VM]
%-----------------------------------
%   Criteria      |  XCell  |  YCell
%-----------------|---------|-------
%            Area |  FEF    |   SC
% [Vis, Mov, Fix] | [1,~,~] | [1,~,~]
%-----------------------------------
%
%      RF of Cells          |
%---------------------------|------------
%         Targ. IN X & IN Y |   AND(X,Y)
%       Targ. IN X NOT in Y |   XOR(X,Y)
%       Targ. IN Y NOT in X |   XOR(X,Y)
% Targ. NOT in X & NOT in Y |   NOT(X|Y)
%-----------------------------------
%    
%%
warning('off');

%% Options for JPSTH computation
binWidth = 1;% use 1 ms for JPSTH computation
% -25 to +25 ms
coincidenceBins = 25;
% area1 = 'SEF';
% area2 = 'FEF';
% area2 = 'SC';

rootAnalysisDir = 'dataProcessed/analysis/JPSTH';
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
for s = 1:size(jpsthCellPairs,1)
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
        if isempty(selTrials)
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
    end % for conditions
    %% Save for the current pair
    tempConditions.cellPairInfo = jpsthPair;
    tempConditions.GithubRef = 'https://github.com/chenchals/thomas-jpsth.git';
    oFn = fullfile(jpsthResultsDir,[pairFilename '.mat']);
    fprintf('Saving processed pair : %s\n',oFn);
    save(oFn,'-v7.3','-struct','tempConditions');
end
end

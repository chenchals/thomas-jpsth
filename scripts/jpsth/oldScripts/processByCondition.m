function processByCondition()

saveFig = true;
%
binWidth = 5;% for JPSTH computations
coincidenceBins = 5;
rootDataDir = '/Volumes/schalllab/data';
rootAnalysisDir = '/Volumes/schalllab/Users/Chenchal/JPSTH';
jpsthResultsDir = fullfile(rootAnalysisDir,['FigsBinWidth' num2str(binWidth,'_%dB')]);
if ~exist(jpsthResultsDir, 'dir')
    mkdir(jpsthResultsDir);
end

% Info files
jpshPairsFile = fullfile(rootAnalysisDir,'JPSTH_PAIRS_CellInfoDB.mat');
trialTypesFile = fullfile(rootAnalysisDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(rootAnalysisDir,'TrialEventTimesDB.mat');
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignEvents = {'CueOn';'SaccadePrimary';'RewardOn'};
% aligned on CueOn
timeWin.(alignEvents{1}) = [-700 400];
% aligned on SaccadePrimary
timeWin.(alignEvents{2}) = [-300 300];
% aligned on RewardOn
timeWin.(alignEvents{3})= [-400 400];
%% Load all JPSTH pair information
% load variable: JpsthPairsCellInfo
jpsthCellPairs = load(jpshPairsFile);
jpsthCellPairs = jpsthCellPairs.JpsthPairCellInfoDB;
jpsthCellPairs.folder = cellfun(@(x) fullfile(rootDataDir,...
    regexprep(x(1),{'D','E'},{'Darwin','Euler'}),'SAT/Matlab'),...
    jpsthCellPairs.datafile,'UniformOutput',false);
sessFiles = jpsthCellPairs.datafile;
rowIdsOfPairsBySession = arrayfun(@(x) find(contains(sessFiles,x)), unique(sessFiles),'UniformOutput',false);
sessions = regexprep(unique(sessFiles),'-RH_.*mat','');

%% Process each cell pair for all available conditions
% for each condition show JPSTH after aligning on all timeWIns specified
% TrialTypes
trialTypes = load(trialTypesFile);
trialTypes = trialTypes.TrialTypesDB;
%retain only sesisons which have jpsth pairs
trialTypes = trialTypes(cellfun(@(x) find(strcmpi(trialTypes.session, x)),sessions),:);
% (Accurate|Fast)*(Correct|ErrorHold|ErrorChoice|ErrorTiming|ErrorNoSaccade)
availConditions = regexp(trialTypes.Properties.VariableNames,'(Accurate.+)|(Fast.+)','match');
availConditions = [availConditions{:}]';
% TrialEventTimes
trialEventTimes = load(trialEventTimesFile);
trialEventTimes = trialEventTimes.TrialEventTimesDB;
%retain only sesisons which have jpsth pairs
trialEventTimes = trialEventTimes(cellfun(@(x) find(strcmpi(trialEventTimes.session, x)),sessions),:);
% (CueOn|FixAcquisition|TargetDeadline|SaccadePrimaryTempo|ToneOn|RewardOn|SaccadePrimary)
% availEventTimes = regexp(trialEventTimes.Properties.VariableNames,'.*$(?<!^session)','match');

%% process all paris by session
for s = numel(rowIdsOfPairsBySession):-1:1
    pairsTodo = jpsthCellPairs(rowIdsOfPairsBySession{s},:);
    nPairs = size(pairsTodo,1);
    file2load = fullfile(pairsTodo.folder{1},pairsTodo.datafile{1});
    [~,sessionName] = fileparts(file2load);
    sessionTrialEventTimes = trialEventTimes(s,:);
    sessionTrialTypes = trialTypes(s,:);
    fprintf('\nDoing JPSTH for session [%s].......\n',sessionName);
    % for each pair of cells
    tempConditions = struct();
    for pair = nPairs:-1:1
        currPair = pairsTodo(pair,:);
        XCellId = currPair.X_cellIdInFile{1};
        YCellId = currPair.Y_cellIdInFile{1};
        pairFilename = char(join({currPair.Pair_UID{1},sessionName,...
            XCellId,currPair.X_area{1},...
            YCellId,currPair.Y_area{1}},'_'));
        fprintf('Processing Pair : %s...\n',pairFilename);
        units = load(file2load,XCellId,YCellId);
        for cond = 1:numel(availConditions)
            condition = availConditions{cond};
            trialNosByCondition = find(sessionTrialTypes.(condition){1});
            tempJpsth = table();
            if isempty(trialNosByCondition)
                tempConditions.(condition) = [];
            else
                for eventId = 1:numel(alignEvents)
                    alignedEvent = alignEvents{eventId};
                    alignedTimeWin = timeWin.(alignedEvent);
                    alignTime = sessionTrialEventTimes.CueOn{1};
                    if isempty(strcmpi(alignedEvent,'CueOn'))
                        alignTime = alignTime + sessionTrialEventTimes.(alignedEvent){1};
                    end
                    XAligned = SpikeUtils.alignSpikeTimes(units.(XCellId),alignTime, alignedTimeWin);
                    YAligned = SpikeUtils.alignSpikeTimes(units.(YCellId),alignTime, alignedTimeWin);
                    
                    XAligned = XAligned(trialNosByCondition,:);
                    YAligned = YAligned(trialNosByCondition,:);
                    
                    temp = SpikeUtils.jpsth(XAligned, YAligned, alignedTimeWin, binWidth, coincidenceBins);
                    tempJpsth(eventId,:) = struct2table(temp,'AsArray',true);
                    %jer = SpikeUtils.jeromiahJpsth(XAligned, YAligned, alignedTimeWin, binWidth, coincidenceBins);                    
                    opts(eventId,1).xCellSpikeTimes = {XAligned}; %#ok<*AGROW>
                    opts(eventId,1).yCellSpikeTimes = {YAligned};
                    opts(eventId,1).trialNosByCondition = {trialNosByCondition};
                    opts(eventId,1).alignedEvent = {alignedEvent};
                    opts(eventId,1).alignedTimeWin = {alignedTimeWin};
                    opts(eventId,1).alignTime = {alignTime(trialNosByCondition)};
                    opts(eventId,1).binWidth = binWidth;
                    opts(eventId,1).coincidenceBins = coincidenceBins;
                end
                
                tempJpsth.Properties.RowNames = alignEvents;
                tempConditions.(condition) = [tempJpsth struct2table(opts,'AsArray',true)];
            end
        end % for condition
        % Save and do next Pair
        tempConditions.cellPairInfo = currPair;
        oFn = fullfile(jpsthResultsDir,[pairFilename '.mat']);
        fprintf('Saving processed pair : %s\n',oFn);
        save(oFn,'-v7.3','-struct','tempConditions');
    end
end %for session group
end
function jpsthProcessByCondition()

    saveFig = true;    
    rootDataDir = '/Volumes/schalllab/data';
    rootAnalysisDir = '/Volumes/schalllab/Users/Chenchal/JPSTH';
    jpsthResultsDir = fullfile(rootAnalysisDir,'FigsBinWidth_1');
    if ~exist(jpsthResultsDir, 'dir')
        mkdir(jpsthResultsDir);
    end
    %
    binWidth = 1;% for JPSTH computations
    coincidenceBinWidth = 5;
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
    pairUids = jpsthCellPairs.Pair_UID;
    rowIdsOfPairsBySession = arrayfun(@(x) find(contains(sessFiles,x)), unique(sessFiles),'UniformOutput',false);
    sessions = regexprep(unique(sessFiles),'-RH_.*mat','');

    %% Process each cell pair for all available conditions
    % for each condition show JPSTH after aligning on all timeWIns specified
    % TrialTypes
    trialTypes = load(trialTypesFile);
    trialTypes = trialTypes.TrialTypesDB;
    % retain only sesisons which have jpsth pairs
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
    availEventTimes = regexp(trialEventTimes.Properties.VariableNames,'.*$(?<!^session)','match');
    availEventTimes = [availEventTimes{:}]';
    % Other variables to load from datafile
    loadOtherVars = [];

    %% process all paris by session
    for s = 1:numel(rowIdsOfPairsBySession)
        sessionJpsths = struct();
        pairsTodo = jpsthCellPairs(rowIdsOfPairsBySession{s},:);   
        unitIdsInFile = unique([pairsTodo.X_cellIdInFile;pairsTodo.Y_cellIdInFile]);
        file2load = fullfile(pairsTodo.folder{1},pairsTodo.datafile{1});
        [~,sessionName] = fileparts(file2load);
        fprintf('\nDoing JPSTH for session [%s].......\n',sessionName);
        vars2load = [loadOtherVars;unitIdsInFile];
        S = load(file2load,vars2load{:});   
        nTrials = size(trialEventTimes.CueOn{s},1);
        nUnits = numel(unitIdsInFile);  
        % Align all spike times Event times to the align Events  
        for a = 1:numel(alignEvents)
            alignEvent = alignEvents{a};
            alignedTimeWin = timeWin.(alignEvent); 
            if strcmpi(alignEvent,'CueOn')
                alignTime = trialEventTimes.CueOn{s};
            else
                alignTime = trialEventTimes.CueOn{s} + trialEventTimes.(alignEvent){s};
            end       
            alignedSpkTimes = cell(nTrials,nUnits);
            for u = 1:nUnits
                unit = S.(unitIdsInFile{u});
                alignedSpkTimes(:,u) = arrayfun(@(x) unit(x,unit(x,:)~=0)-alignTime(x),(1:nTrials)','UniformOutput',false);
                alignedSpkTimes(:,u) = cellfun(@(x) x(x>=alignedTimeWin(1) & x<=alignedTimeWin(2)),alignedSpkTimes(:,u),'UniformOutput',false);
                %alignedSpkTimes(:,u) = temp;
            end
            for cond = 1:numel(availConditions)
                condition = availConditions{cond};
                trialNosByCondition = find(trialTypes.(condition){s}); 
                if isempty(trialNosByCondition)
                    sessionJpsths.(condition).(alignEvent) = [];
                    continue;
                end
                fprintf('Doing JPSTH for condition [%s] - aligned on event [%s]...',condition,alignEvent);
                [~,jpsthTable] = newJpsth(alignedSpkTimes(trialNosByCondition,:),unitIdsInFile,alignedTimeWin,binWidth,coincidenceBinWidth);
                %verify cell-pairing in newJpsh with pairsTodo
                jpsthStruct_pairKeys = strcat(jpsthTable.xCellId,'-',jpsthTable.yCellId);
                pairsTodo_pairKeys = strcat(pairsTodo.X_cellIdInFile,'-',pairsTodo.Y_cellIdInFile);
                % the number of pairs match AND all the rows match as-is or
                % one of the arrays conatin all items may be in a different order?
                if isequaln(jpsthStruct_pairKeys,pairsTodo_pairKeys) || ...
                        sum(contains(jpsthStruct_pairKeys,pairsTodo_pairKeys)) == numel(jpsthStruct_pairKeys)
                    jpsthTable.sortIdxs = cellfun(@(x) find(contains(pairsTodo_pairKeys,x)),jpsthStruct_pairKeys);
                    jpsthTable =  [pairsTodo sortrows(jpsthTable,'sortIdxs')];
                else
                    error('*****Number of pairs done by newJpsth does not match the pairsToDo\n******\n');
                end
                jpsthTable.trialNosByCondition(:,1) = {trialNosByCondition};
                sessionJpsths.(condition).(alignEvent) = jpsthTable;            
            end
        end
        if saveFig
            % Save all jpsth pairs for session
            fprintf('Saving JPSTH fo all pairs for session [%s]\n',sessionName); %#ok<UNRCH>
            save(fullfile(jpsthResultsDir,[sessionName '_JPSTH']),'-v7.3','-struct','sessionJpsths');
        else
            fprintf('Done JPSTH fo all pairs for session [%s]\n',sessionName);
        end
    end
    fprintf('Done...\n');
end
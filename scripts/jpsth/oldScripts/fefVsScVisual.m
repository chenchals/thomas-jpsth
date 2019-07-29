%% Cross-area JPSTH analysis of cell pairs
% FEF_SC both X and Y cell must be Visual 
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
%% Options for JPSTH computation
binWidth = 1;% use 1 ms for JPSTH computation
% -25 to +25 ms
coincidenceBins = 25;
area1 = 'FEF';
area2 = 'SC';

rootDataDir = '/Volumes/schalllab/data';
rootAnalysisDir = 'dataProcessed/JPSTH';
jpsthResultsDir = fullfile(rootAnalysisDir,['FEF_SC_Visual' num2str(binWidth,'_%dms')]);
if ~exist(jpsthResultsDir, 'dir')
    mkdir(jpsthResultsDir);
end

% Info files
jpshPairsFile = fullfile(rootAnalysisDir,'JPSTH_PAIRS_CellInfoDB.mat');
trialTypesFile = fullfile(rootAnalysisDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(rootAnalysisDir,'TrialEventTimesDB.mat');
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignEventTimeWin = containers.Map;
alignEventTimeWin('CueOn') = [-200 400];
alignEventTimeWin('SaccadePrimary') = [-300 300];
alignEventTimeWin('RewardOn') = [-300 300];
alignEvents = alignEventTimeWin.keys;
%% Load all JPSTH pair information
% load variable: JpsthPairsCellInfo
jpsthCellPairs = load(jpshPairsFile);
jpsthCellPairs = jpsthCellPairs.JpsthPairCellInfoDB;

% Filter cell pairs for FEF/SC
jpsthCellPairs = jpsthCellPairs(contains(jpsthCellPairs.X_area,area1) & contains(jpsthCellPairs.Y_area,area2),:);
% Filter cell pairs for Visual = 1.0
visIdx = cellfun(@(x,y) x(1)==1 & y(1)==1, jpsthCellPairs.X_visMovFix,jpsthCellPairs.Y_visMovFix);
jpsthCellPairs = jpsthCellPairs(visIdx,:);

% Data files for the filtered cell pairs
jpsthCellPairs.folder = cellfun(@(x) fullfile(rootDataDir,...
    regexprep(x(1),{'D','E'},{'Darwin','Euler'}),'SAT/Matlab'),...
    jpsthCellPairs.datafile,'UniformOutput',false);
sessFiles = jpsthCellPairs.datafile;
% Group by session
rowIdsOfPairsBySession = arrayfun(@(x) find(contains(jpsthCellPairs.datafile,x)), ...
                    unique(jpsthCellPairs.datafile),'UniformOutput',false);
sessions = regexprep(unique(jpsthCellPairs.datafile),'-RH_.*mat','');

%% Load Trial Types and Trial Event Times and filter for the sessions of interest above
% TrialTypes
trialTypes = load(trialTypesFile);
trialTypes = trialTypes.TrialTypesDB;
% Filter trialTypes to have only sesisons of interest
trialTypes = trialTypes(cellfun(@(x) find(strcmpi(trialTypes.session, x)),sessions),:);
% TrialEventTimes
trialEventTimes = load(trialEventTimesFile);
trialEventTimes = trialEventTimes.TrialEventTimesDB;
%retain only sesisons which have jpsth pairs
trialEventTimes = trialEventTimes(cellfun(@(x) find(strcmpi(trialEventTimes.session, x)),sessions),:);
% Available conditions - 
% (Accurate|Fast)*(Correct|ErrorHold|ErrorChoice|ErrorTiming|ErrorNoSaccade)
availConditions = regexp(trialTypes.Properties.VariableNames,'(Accurate.+)|(Fast.+)','match');
availConditions = [availConditions{:}]';

%% Group by singleton location in RF fo X and Y cell
rfLocNames = {'TargetInXandY','TargetInXnotY','TargetInYnotX','TargetNotInXorY'};
fx_groupRFs = @(xRF,yRF) deal(xRF, yRF, intersect(xRF,yRF), setdiff(xRF,yRF),...
                setdiff(yRF,xRF), setdiff(0:7,[xRF(:);yRF(:)]));
[temp.xRF,temp.yRF,temp.(rfLocNames{1}),temp.(rfLocNames{2}),temp.(rfLocNames{3}),temp.(rfLocNames{4})] = ...
    cellfun(@(xRF,yRF) fx_groupRFs(xRF,yRF), jpsthCellPairs.X_RF, jpsthCellPairs.Y_RF,'UniformOutput',false);
rfLocsGroups = struct2table(temp); clear temp;
%% For each grouped cell pairs by session dir JPSTHs for every pair
for s = 1:numel(rowIdsOfPairsBySession)
    rowIdsForPairs = rowIdsOfPairsBySession{s};
    pairsTodo = jpsthCellPairs(rowIdsForPairs,:);
    nPairs = size(pairsTodo,1);
    file2load = fullfile(pairsTodo.folder{1},pairsTodo.datafile{1});
    [~,sessionName] = fileparts(file2load);
    sessionTrialEventTimes = trialEventTimes(s,:);
    sessionTrialTypes = trialTypes(s,:);
    sessionRfLocs = rfLocsGroups(rowIdsForPairs,:);
    fprintf('\nDoing JPSTH for session [%s].......\n',sessionName);
    % for each pair of cells
    for pair = 1:nPairs
        tempConditions = struct();
        currPair = pairsTodo(pair,:);
        XCellId = currPair.X_cellIdInFile{1};
        YCellId = currPair.Y_cellIdInFile{1};
        pairFilename = char(join({currPair.Pair_UID{1},sessionName,...
            XCellId,currPair.X_area{1},...
            YCellId,currPair.Y_area{1}},'_'));
        fprintf('Processing Pair : %s...\n',pairFilename);
        units = load(file2load,XCellId,YCellId);
        % for each named singleton location
        for rf = 1:numel(rfLocNames)
            rfLocName = rfLocNames{rf};
            rfLocs = sessionRfLocs.(rfLocName){pair};
            if isempty(rfLocs) || any(isnan(rfLocs))
                tempConditions.(rfLocName) = [];
            else
                % for each condition (outcome)
                for cond = 1:numel(availConditions)                    
                    condition = availConditions{cond};
                    % Get trial Nos for rfLoc and condition
                    selTrials = find(sessionTrialTypes.(condition){:} & ismember(sessionTrialTypes.SingletonLoc{:},rfLocs));                    
                    if isempty(selTrials)
                        tempConditions.(rfLocName).(condition) = [];
                        continue;
                    end
                    % for each aligned event
                    for evId = 1:numel(alignEvents)
                        alignedEvent = alignEvents{evId};
                        alignedTimeWin = alignEventTimeWin(alignedEvent);
                        alignTime = sessionTrialEventTimes.CueOn{1};
                        if isempty(strcmpi(alignedEvent,'CueOn'))
                            alignTime = alignTime + sessionTrialEventTimes.(alignedEvent){1};
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
                        opts(evId,1).alignedEvent = {alignedEvent};
                        opts(evId,1).alignedTimeWin = {alignedTimeWin};
                        opts(evId,1).alignTime = {alignTime};
                        opts(evId,1).binWidth = binWidth;
                        opts(evId,1).coincidenceBins = coincidenceBins;                       
                    end % for alignEvents
                    tempJpsth.Properties.RowNames = alignEvents;
                    tempConditions.(rfLocName).(condition) = [tempJpsth struct2table(opts,'AsArray',true)];                    
                end % for conditions
            end   
        end % for rfLocNames
        % Save for the current pair
        tempConditions.cellPairInfo = currPair;
        tempConditions.singletonLocs = sessionRfLocs(pair,:);
        tempConditions.GithubRef = 'git clone --branch JPSTH_FEFxSC_V1 https://github.com/chenchals/schalllab-jpsth.git';
        oFn = fullfile(jpsthResultsDir,[pairFilename '.mat']);
        fprintf('Saving processed pair : %s\n',oFn);
        save(oFn,'-v7.3','-struct','tempConditions');
    end
end



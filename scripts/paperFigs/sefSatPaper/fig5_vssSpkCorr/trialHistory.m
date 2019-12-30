% Create table for outcomes by session split on previous outcome.
% Previous outcome criteria:
%   1. Previous outcome must be within the same block as current outcome
%   2. Previous outcome must not be the first oucome in the block
% Criteria(1) ensures that the monk is in a given state of either Fast or
%      Accurate as these SAT conditions are blocked 
% Criteria(2) ensures that the monk is not "transitioning" from one state
% to another.
% The possible combinations are the foillowing:
% cOC: current Outcome
% pOC: previous outcome [1,2 or 3]
% A : Accurate,  F: Fast
% _C : Correct, _EC: ErrorChoice, _ET: ErrorTiming
% _____________________________
% | cOC  | pOC1 | pOC2 | pOC3 |  
% |------|------|------|------|
% | A_C  | A_C  | A_EC | A_ET |
% | A_EC | A_C  | A_EC | A_ET |
% | A_ET | A_C  | A_EC | A_ET |
% |------|------|------|------|
% | F_C  | F_C  | F_EC | F_ET |
% | F_EC | F_C  | F_EC | F_ET |
% | F_ET | F_C  | F_EC | F_ET |
% -----------------------------
%

% load the behavior file for getting trial outcomes
trialHistoryFile = 'dataProcessed/analysis/spkCorr/trialHistory/SpkCorr_TrialHistory.mat';
trialTypeDbFile = 'dataProcessed/dataset/TrialTypesDB.mat';
trialTypesDB = load(trialTypeDbFile);
trialTypesDB = trialTypesDB.TrialTypesDB;

trialTypesDB.monk = cellfun(@(x) x(1),trialTypesDB.session);
trialTypesDB = trialTypesDB(ismember(trialTypesDB.monk,{'D','E'}),:);

% colNames =  trialTypesDB.Properties.VariableNames';
%
% fastCols = colNames(contains(colNames,'Fast'));
% accCols = colNames(contains(colNames,'Accurate'));
%% Use trials as in beh - ignore block boundaries 
fastCols = {'FastCorrect','FastErrorChoice','FastErrorTiming'}';
accCols = {'AccurateCorrect','AccurateErrorChoice','AccurateErrorTiming'}';
useCols = [fastCols;accCols];
%%
nthPrevTrial = 1;
out = struct();
[out.(['trialHistory_' num2str(nthPrevTrial,'%d')])] = getNthTrialHistoryUseBlock(trialTypesDB,useCols,nthPrevTrial);

if exist(trialHistoryFile,'file')
    save(trialHistoryFile,'-append','-v7.3','-struct','out');
else
    save(trialHistoryFile,'-v7.3','-struct','out');
end


%%
function [prevHistoryTbl] = getNthTrialHistoryUseBlock(trialTypesDB,useColsNames,nthPrevTrial)
nPrevTrial = nthPrevTrial; % number of trials before current trial
prevHistoryTbl = table();
useColsNames = useColsNames';
for s = 1:size(trialTypesDB.session)
    tbl = table();
    sess = trialTypesDB.session{s};
    fprintf('Doing Session : %s\n',sess);
    tbl.TrialBlockNum = trialTypesDB.TrialBlockNum{1};
    for uc = useColsNames
       tbl.(uc{1}) = trialTypesDB.(uc{1}){1};
    end
    % Where does SAT trial start
    tbl.TrialBlockNum(isnan(tbl.TrialBlockNum)) = 0;
    % satOffsetIdx = find(isnan(tbl.TrialBlockNum), 1, 'last' );
    %tbl = tbl(satOffsetIdx + 1:end,:);
    % max number of blocks in session
    maxBlkNo = nanmax(tbl.TrialBlockNum);
    % within each block accumulate outcome of prev trial for current trl
    % outcome 
    tempBlk = table(); 
    parfor blk = 1:maxBlkNo
        blkStartIdx = find(tbl.TrialBlockNum==blk, 1, 'first' );
        blkEndIdx = find(tbl.TrialBlockNum==blk, 1, 'last' );
        % do not use 1st trial in the block due to STATE change from FAST
        % to ACCURATE or vice versa
        blkStartIdx = blkStartIdx + 1;
        % nth trial in block to last trial in blk
        currTrlsIdx = (blkStartIdx + nPrevTrial:blkEndIdx);
        prevTrlsIdx = (blkStartIdx:blkEndIdx - nPrevTrial);
        tempOutcome = table();
        for oCome = 1:numel(useColsNames)
            currOutcome = useColsNames{oCome};
            prevOutcomes = useColsNames;            
            prevOutcomesColNames = strcat(prevOutcomes,num2str(nPrevTrial,'_p%d'));  
            prevOutcomesCountColNames = strcat(prevOutcomes,num2str(nPrevTrial,'_p%d_count'));  
            temp = table();          
            temp.session = {sess};
            temp.blkNum = blk;
            % trls with cirrOutcome in the block
            trlsCurrOutcome = tbl.(currOutcome)(currTrlsIdx) == 1; 
            if ~any(trlsCurrOutcome)
                % continue ?
                temp.currOutcome = {currOutcome};
                temp.currTrials = {[]};
                temp.currTrialsCount = 0;
                temp{1,prevOutcomesColNames} = repmat({[]},1,numel(prevOutcomesColNames));
                temp(1,prevOutcomesCountColNames) = repmat({0},1,numel(prevOutcomesCountColNames));
                tempOutcome = [tempOutcome; temp];
                continue;
            end           
            trlsPrevGivenCurrOutcome = prevTrlsIdx(trlsCurrOutcome);
            
            % Trial nos of prev outcome, given current outcome
            trlsPrevOutcomeGivenCurrOutcome = cellfun(@(x) {trlsPrevGivenCurrOutcome(tbl.(x)(trlsPrevGivenCurrOutcome) == 1)}, prevOutcomes,'UniformOutput',false);
            trlsPrevOutcomeGivenCurrOutcome = cell2table(trlsPrevOutcomeGivenCurrOutcome,'VariableNames',prevOutcomesColNames);
            % Trial counts of prev outcome given current outcome
            trlsCountPrevOutcomeGivenCurrOutcome = cellfun(@(x) numel(trlsPrevGivenCurrOutcome(tbl.(x)(trlsPrevGivenCurrOutcome) == 1)), prevOutcomes,'UniformOutput',false);
            trlsCountPrevOutcomeGivenCurrOutcome = cell2table(trlsCountPrevOutcomeGivenCurrOutcome,'VariableNames',prevOutcomesCountColNames);
            % aggregate

            temp.currOutcome = {currOutcome};
            temp.currTrials = {currTrlsIdx(trlsCurrOutcome)};
            temp.currTrialsCount= numel(currTrlsIdx(trlsCurrOutcome));
            tempOutcome = [tempOutcome; [temp trlsPrevOutcomeGivenCurrOutcome trlsCountPrevOutcomeGivenCurrOutcome]];
            
        end
        tempBlk = [tempBlk;tempOutcome];
    end
    sessTbl = aggregateBlocks(tempBlk);
    prevHistoryTbl = [prevHistoryTbl;sessTbl]; %#ok<*AGROW>
end

end

function [sessTbl] = aggregateBlocks(tempBlk)
    % inner function
    warning('off')
    colNames = tempBlk.Properties.VariableNames;
    prevTrlsCols = colNames(~cellfun(@isempty,regexp(colNames,'.*_p\d*$','match')));
    prevTrlsCountCols = colNames(~cellfun(@isempty,regexp(colNames,'.*_count\d*$','match')));

    currOutcomes = unique(tempBlk.currOutcome,'stable');
    sessTbl = table();
    for co = 1:numel(currOutcomes)
        currOutcome = currOutcomes{co};
        idx = find(ismember(tempBlk.currOutcome,currOutcome));
        sessTbl.session{co} = tempBlk.session{idx(1)};
        sessTbl.currOutcome{co} = currOutcome;
        currTrls = [tempBlk.currTrials{idx}];
        sessTbl.currOutcomeTrls{co} = currTrls;
        sessTbl.currOutcomeTrlsCount(co) = sum(tempBlk.currTrialsCount(idx));
        % trial counts
        sessTbl(co,prevTrlsCountCols) = cellfun(@(x) sum(tempBlk.(x)(idx)),prevTrlsCountCols,'UniformOutput',false);
        % trialNos
        sessTbl{co,prevTrlsCols} = cellfun(@(x) [tempBlk.(x){idx}],prevTrlsCols,'UniformOutput',false);
    end
end


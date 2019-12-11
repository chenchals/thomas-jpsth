% Extract end points of all saccades from movesAll data
behFn = 'dataProcessed/dataset/dataBehavior_SAT_DaEu.mat';
movesAllFn = 'dataProcessed/dataset/movesAll_SAT.mat';

behDaEuIdx ={ 1:8;9:14};
beh = load(behFn);
binfoDaEu = beh.binfoSAT_DaEu;

saccades_Da = table();
saccades_Eu = table();

movesAll = load(movesAllFn);

moveDaEu = {'movesAllDa','movesAllEu'};
moveFields = {'index' 'octant' 'resptime' 'trial' 'duration'};
warning('off')
outTbl = table();
tic
for mm = 1:2
    sessIdxs = behDaEuIdx{mm};
    binfo = binfoDaEu(sessIdxs,:);
    moves = struct2table(movesAll.(moveDaEu{mm}));
    parfor s = 1:numel(sessIdxs)
        tic
        temp = table();
        currBeh = binfo(s,:);
        currMoves = moves(s,:);
        behTrialNums = (1:currBeh.num_trials)';       
        % extract for session
        saccTrialNums = double(sort(unique(currMoves.trial{1})))';
        missingTrialNums = setdiff(behTrialNums,saccTrialNums);
        
        temp.trialNum = saccTrialNums;
        temp.nSaccades = arrayfun(@(x) sum(currMoves.trial{1}==x & currMoves.index{1}>0), saccTrialNums);
        temp.saccIndex = arrayfun(@(x) currMoves.index{1}(currMoves.trial{1}==x), saccTrialNums,'UniformOutput',false);
        temp.resptime = arrayfun(@(x) double(currMoves.resptime{1}(currMoves.trial{1}==x)), saccTrialNums,'UniformOutput',false);
        temp.saccDuration = arrayfun(@(x) double(currMoves.duration{1}(currMoves.trial{1}==x)), saccTrialNums,'UniformOutput',false);
        temp.saccOctant = arrayfun(@(x) double(currMoves.octant{1}(currMoves.trial{1}==x)), saccTrialNums,'UniformOutput',false);
        
        % add missing trials and sort on trial number
        if ~isempty(missingTrialNums)
            temp.trialNum(end+1:end+numel(missingTrialNums)) = missingTrialNums;
        end
        temp = sortrows(temp,{'trialNum'});        
        % add trgat position        
        temp.targetOctant = double(currBeh.tgt_octant{1}');
        
        % saccIdxOnTarget - which saccade landed on target 2? 3rd? ...
        temp.saccIndexOnTarget = arrayfun(@(x) max([find(temp.saccOctant{x}==temp.targetOctant(x)),NaN]),temp.trialNum);
          
        % create out table
        fns = temp.Properties.VariableNames;
        ot = table();
        for f = 1:numel(fns)
            ot.(fns{f}) = {temp.(fns{f})};
        end
        ot.session = currBeh.Properties.RowNames;
        ot.num_trials = currBeh.num_trials;
        ot.Properties.RowNames = currBeh.Properties.RowNames;
        % aggregate all sessions
        outTbl = [outTbl;ot]; 
        fprintf('Done session %s (%4.2f)\n',currBeh.Properties.RowNames{1},toc)
    end
end
toc



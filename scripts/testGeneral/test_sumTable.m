tempBlk = load('tempBlk.mat');
tempBlk = tempBlk.tempBlk;
colNames = tempBlk.Properties.VariableNames;
prevTrlsCols = colNames(~cellfun(@isempty,regexp(colNames,'.*_p\d*$','match')));
prevTrlsCountCols = colNames(~cellfun(@isempty,regexp(colNames,'.*_count\d*$','match')));

currOutcomes = unique(tempBlk.currOutcome,'stable');
sessTbl = table();
for co = 1:numel(currOutcomes)
    currOutcome = currOutcomes{co};
    idx = find(ismember(tempBlk.currOutcome,currOutcome));
    currTbl = tempBlk(idx,:);
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
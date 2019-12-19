% load the behavior file for getting trial outcomes

trialTypeDbFile = 'dataProcessed/dataset/TrialTypesDB.mat';
trialTypesDB = load(trialTypeDbFile);
trialTypesDB = trialTypesDB.TrialTypesDB;

trialTypesDB.monk = cellfun(@(x) x(1),trialTypesDB.session);
trialTypesDB = trialTypesDB(ismember(trialTypesDB.monk,{'D','E'}),:);

colNames =  trialTypesDB.Properties.VariableNames';

fastCols = colNames(contains(colNames,'Fast'));
accCols = colNames(contains(colNames,'Accurate'));

%%
prevFastAcc = table();
prevFastAccCounts = table();
useCols = [fastCols;accCols]';
for s = 1:size(trialTypesDB.session)
    tbl = trialTypesDB(s,:);
    sess = tbl.session{1};
    for ii = 1:numel(useCols)
        temp = table();
        outcome_n = useCols{ii};
        prevIdx = circshift(tbl.(outcome_n){1},-1) == 1;
        prevTrls = cellfun(@(x) find(prevIdx & (tbl.(x){1} == 1)),useCols,'UniformOutput',false);
        prevTrls = cell2table(prevTrls,'VariableNames',useCols);
        prevTrlsCounts = cellfun(@(x) numel(find(prevIdx & (tbl.(x){1} == 1))),useCols,'UniformOutput',false);
        prevTrlsCounts = cell2table(prevTrlsCounts,'VariableNames',useCols);
        temp.session = {sess};
        temp.currOutcome = {outcome_n};
        prevFastAcc = [prevFastAcc;[temp prevTrls]];
        prevFastAccCounts = [prevFastAccCounts;[temp prevTrlsCounts]];
    end
end
%%
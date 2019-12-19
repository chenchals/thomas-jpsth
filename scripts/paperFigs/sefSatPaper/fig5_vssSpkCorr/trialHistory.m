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
nPrevTrial = 1; % number of trial before current trial
prevColNames = strcat(useCols,num2str(nPrevTrial,'_p%d'));
for s = 1:size(trialTypesDB.session)
    tbl = trialTypesDB(s,:);
    satStart = find(isnan(tbl.Fast{1}), 1, 'last' );
    sess = tbl.session{1};
    for ii = 1:numel(useCols)
        temp = table();
        outcome_n = useCols{ii};
        % from trial#2 to end
        idxOffset = satStart + nPrevTrial + 1;        
        prevIdx = tbl.(outcome_n){1}(idxOffset+1:end) == 1;
        prevTrls = cellfun(@(x) find(prevIdx & (tbl.(x){1}(1:end-idxOffset) == 1)),useCols,'UniformOutput',false);
        prevTrls = cell2table(prevTrls,'VariableNames',prevColNames);
        prevTrlsCounts = cellfun(@(x) numel(find(prevIdx & (tbl.(x){1}(1:end-idxOffset) == 1))),useCols,'UniformOutput',false);
        prevTrlsCounts = cell2table(prevTrlsCounts,'VariableNames',prevColNames);
        temp.session = {sess};
        temp.currOutcome = {outcome_n};
        prevFastAcc = [prevFastAcc;[temp prevTrls]];
        prevFastAccCounts = [prevFastAccCounts;[temp prevTrlsCounts]];
    end
end
%%
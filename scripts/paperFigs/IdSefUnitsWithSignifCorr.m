% Identify SEF units that have significant correlation with more than 1 FEF
% and or SC units

%% Load static RSC matrix for all pairs of SEF units
    fnUnitInfo = 'dataProcessed/dataset/dataNeurophys_SAT.mat';% for unit info
    fnRho = 'dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticRhoPval.mat';
        % Load data
    temp = load(fnRho,'-regexp'); % load cross areas and self areas
    temp.spkCorrColumnDefs=[];
    allSpkCorr = table();
    fns = fieldnames(temp);
    for jj = 1:numel(fns)
        allSpkCorr = [allSpkCorr;temp.(fns{jj})]; %#ok<*AGROW>
    end
    unitInfo = load(fnUnitInfo,'unitInfo');
    unitInfo = unitInfo.unitInfo;
    
%% Summarize the spike corr table for PostSaccade by outcome by number of signifcant connections for a SEF unit
% Filter spkCorr for only SEF-FEF, SEF-SC connections
spkCorr = allSpkCorr(ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,{'FEF','SC'}),:);
% Filter for PostSaccade
epoch = 'PostSaccade';
spkCorr = spkCorr(ismember(spkCorr.alignedName,epoch),:);
% filter for significant Rsc
spkCorr = spkCorr(spkCorr.pvalRaw_150ms <= 0.01,:);
% get unique SEF unit list
sefUnits = unique(spkCorr.X_unitNum);
satConds = {'Fast','Accurate'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
outTable = table();
for u = 1:numel(sefUnits)
    unitNum = sefUnits(u);
    uInfo = unitInfo(unitInfo.unitNum == unitNum,:);
    for c = 1:numel(satConds)
        satCond = satConds{c};
        for o = 1:numel(outcomes)
            outcome = outcomes{o};
            temp = spkCorr(spkCorr.X_unitNum == unitNum & ...
                           contains(spkCorr.condition,satCond) & ...
                           contains(spkCorr.condition,outcome),:);
            if size(temp,1)<2
                continue;
            end
            t = table();
            
            t.monk{1} = uInfo.monkey;
            t.sessNum(1) = uInfo.sessNum;
            t.sess{1} = uInfo.sess;
            t.unit{1} = uInfo.unit;
            
            t.satCondition{1} = satCond;            
            t.outcome{1} = outcome;
            t.epoch{1} = epoch;

            t.sefUnit(1) = unitNum;
            t.numSignifConn(1) = size(temp,1);
            
            fef = temp(strcmp(temp.Y_area,'FEF'),:);
            sc = temp(strcmp(temp.Y_area,'SC'),:);
            
            t.numFefConn(1) = size(fef,1);
            t.numScConn(1) = size(sc,1);
            
            cu ='';
            if ~isempty(fef.Y_unitNum)
                cu = sprintf('%d,',fef.Y_unitNum);
            end
            cu = ['[',cu(1:end-1),']'];
            t.fefUnits{1} = cu;
            
            cu = '';
            if ~isempty(sc.Y_unitNum)
                cu = sprintf('%d,',sc.Y_unitNum);
            end
            cu = ['[',cu(1:end-1),']'];
            t.scUnits{1} = cu;
            
            outTable = [outTable;t];
        end
    end
end
writetable(outTable,'dataProcessed/analysis/spkCorr/SEFunitsWithGT1SignifRsc.csv')


% Run anova for different factors of SAT for spike correlations
fn = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat';
spkCorr = load(fn);
outFn = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAnovaSummary.mat';
useAbsRho = 1; % no support for useAbsRho = 0
%% aggregate required fields into a single table
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
% Prune to columns of interest for each pairArea
useCols = {
    'pairAreas'
    'XY_Dist'
    'condition'
    'alignedName'
    'rhoRaw_200ms'
    };

spkCorrAllTbl = table();
for pa =1:numel(pairAreas)
    pairArea = pairAreas{pa};
    spkCorrAllTbl = [spkCorrAllTbl; spkCorr.(pairArea)(:,useCols)];
end
%% Filter data and assemble values and factors into a table
% Filter data --> remove all values where XY_Dist = 0
spkCorrAllTbl = spkCorrAllTbl(spkCorrAllTbl.XY_Dist~=0,:);
if useAbsRho
    spkCorrAllTbl.rhoRaw_200ms = abs(spkCorrAllTbl.rhoRaw_200ms);
end

idx_TrialOutcome = ismember(spkCorrAllTbl.condition, {'AccurateErrorChoice','FastErrorChoice'});
spkCorrAllTbl = spkCorrAllTbl(idx_TrialOutcome,:);
%% Recode groups/factors for anova - CONDITION (2) by EPOCH (4) = (8*7)/2 = 28 comparisions
valsGroupsTbl = table();
valsGroupsTbl.yVals = spkCorrAllTbl.rhoRaw_200ms;
valsGroupsTbl.condition = regexprep(spkCorrAllTbl.condition,'Correct|Error.*','');
valsGroupsTbl.epoch = spkCorrAllTbl.alignedName;

[conditionByEpoch] = satAnova(valsGroupsTbl);

%% Recode groups/factors for anova - PAIRAREA (3) by EPOCH (4) = (12*11)/2 = 66 comparisions
valsGroupsTbl = table();
valsGroupsTbl.yVals = spkCorrAllTbl.rhoRaw_200ms;
valsGroupsTbl.pairAreas = spkCorrAllTbl.pairAreas;
valsGroupsTbl.epoch = spkCorrAllTbl.alignedName;

[pairAreasByEpoch] = satAnova(valsGroupsTbl);

%% Recode groups/factors for anova - CONDITION (2) by PAIRAREA (3) by EPOCH (4) = (24*23/2) = 276 comparisions
% valsGroupsTbl = table();
% valsGroupsTbl.yVals = spkCorrAllTbl.rhoRaw_200ms;
% valsGroupsTbl.condition = regexprep(spkCorrAllTbl.condition,'Correct|Error.*','');
% valsGroupsTbl.pairAreas = spkCorrAllTbl.pairAreas;
% valsGroupsTbl.epoch = spkCorrAllTbl.alignedName;
% 
% [conditionByPairAreasByEpoch] = satAnova(valsGroupsTbl);
% 
% 


%%
%save(outFn,'conditionByEpoch','pairAreasByEpoch');





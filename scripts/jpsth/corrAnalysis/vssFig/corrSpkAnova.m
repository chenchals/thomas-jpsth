%corrSpkAnova.m
% This script collects static rSC values for four trial epochs: Baseline,
% Visual Response, Post-Saccade, and Post-Reward. It computes an analysis
% of variance with factors Trial Epoch (4 levels) and Task Condition (2
% levels).

areaTest = 'SEF-FEF';
trialOutcome = {'Correct','ErrorChoice','ErrorTiming'};
% outFn = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAnovaSummary.mat';

%% Run anova for different factors of SAT for spike correlations
spkCorr = load('dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat');

%% Aggregate required fields into a single table
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
% Prune to columns of interest for each pairArea
useCols = {'pairAreas','XY_Dist','condition','alignedName','rhoRaw_200ms'};

spkCorr_All = table();
for pa =1:numel(pairAreas)
    pairArea = pairAreas{pa};
    spkCorr_All = [spkCorr_All; spkCorr.(pairArea)(:,useCols)];
end

%remove all values where XY_Dist = 0 -- commented by TR 2019-12-02
% spkCorr_All = spkCorr_All(spkCorr_All.XY_Dist~=0,:);

%take absolute value of correlation
spkCorr_All.rhoRaw_200ms = abs(spkCorr_All.rhoRaw_200ms);


for to = 1:3

if (USE_UNSIGNED_RHO)
    
end

idxTrialOutcome = ismember(spkCorr_All.condition, {'AccurateErrorChoice','FastErrorChoice'});
idxArea = strcmp(spkCorr_All.pairAreas,areaTest);
spkCorrFilt = spkCorr_All(idxTrialOutcome & idxArea,:);

%% Recode groups/factors for anova - CONDITION (2) by EPOCH (4) = (8*7)/2 = 28 comparisions
valsGroupsTbl = table();
valsGroupsTbl.yVals = spkCorrFilt.rhoRaw_200ms;
valsGroupsTbl.condition = regexprep(spkCorrFilt.condition,'Correct|Error.*','');
valsGroupsTbl.epoch = spkCorrFilt.alignedName;
[conditionByEpoch] = satAnova(valsGroupsTbl);
end % for : trialOutcome (to)

%% Save data
%save(outFn,'conditionByEpoch','pairAreasByEpoch');



% %% Recode groups/factors for anova - PAIRAREA (3) by EPOCH (4) = (12*11)/2 = 66 comparisions
% valsGroupsTbl = table();
% valsGroupsTbl.yVals = spkCorrFilt.rhoRaw_200ms;
% valsGroupsTbl.pairAreas = spkCorrFilt.pairAreas;
% valsGroupsTbl.epoch = spkCorrFilt.alignedName;
% 
% [pairAreasByEpoch] = satAnova(valsGroupsTbl);






%load r_sc data
load('dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticNew.mat', 'SEF_FEF','SEF_SC','SEF_SEF');

%% Prepare to extract data for ANOVA
rhoColName = 'rhoRaw_150ms';
useCols = {'pairAreas','XY_Dist','condition','alignedName',rhoColName};

%isolate area pairing of interest
spkCorr = [SEF_FEF(:,useCols); SEF_SC(:,useCols)];
areaPair = 'SEF-FEF & SEF-SC';

%take absolute value of correlation
spkCorr.rho = abs(spkCorr.(rhoColName));
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.epoch = spkCorr.alignedName;
epochs = {'Baseline','Visual','PostSaccade','PostReward'};

%% Do anova for the following
% The Y value for all Anovas come foprm 
% y = column 'rho' which could be 
%      from: {'rhoRaw_50ms' or 'rhoRaw_150ms' or'rhoRaw_200ms'}
% Anova (1) Hold Factor: 'epoch' constant which comes from column 'alignedName', and do a :
%      2-way anova Factors:
%      F1. 'satCondition' from column 'condition', values: Fast|Accurate
%      F2. 'outcome' from column 'condition', values: Correct|ErrorChoice|ErrorTiming  
statsAnova = struct();
anovaModelName = 'interaction'; %  [linear|interaction]
doMultiCompareFlag = 1; % [0|1] 
alpha = 0.05;
for ep = 1:numel(epochs)
    epoch = epochs{ep};
    idx = ismember(spkCorr.epoch,epoch);
    anovaTbl = table();
    anovaTbl.yVals = spkCorr{idx,'rho'};
    anovaTbl.satCondition = spkCorr{idx,'satCondition'};
    anovaTbl.outcome = spkCorr{idx,'outcome'};
    statsAnova.(epoch) = satAnova(anovaTbl,'interaction',doMultiCompareFlag,alpha);
end



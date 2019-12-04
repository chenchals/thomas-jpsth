%corrSpkAnova.m
% This script collects static rSC values for four trial epochs: Baseline,
% Visual Response, Post-Saccade, and Post-Reward. It computes an analysis
% of variance with factors Trial Epoch (4 levels) and Task Condition (2
% levels).

windowTest = 'rhoRaw_150ms';
trialOutcome = {'Correct','ErrorChoice','ErrorTiming'};
ylimPlot = [0.04 0.18];

%load r_sc data
% load('dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat', 'SEF_FEF','SEF_SC','SEF_SEF');

%% Aggregate required fields into a single table
% pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
% Prune to columns of interest for each pairArea
useCols = {'pairAreas','XY_Dist','condition','alignedName',windowTest};

%isolate area pairing of interest
% spkCorr = SEF_SEF(:,useCols);
spkCorr = [SEF_FEF(:,useCols); SEF_SC(:,useCols)];
areaPair = 'SEF-FEF & SEF-SC';

%take absolute value of correlation
spkCorr.(windowTest) = abs(spkCorr.(windowTest));

%remove all values where XY_Dist = 0 -- commented by TR 2019-12-02
% spkCorr_All = spkCorr_All(spkCorr_All.XY_Dist~=0,:);

%initialize output for plot and ANOVA
rsc_Acc = cell(1,length(trialOutcome));
rsc_Fast = cell(1,length(trialOutcome));
stats_Anova = cell(1,length(trialOutcome));

rsc_Baseline_Acc = NaN(2,length(trialOutcome));
rsc_Baseline_Fast = NaN(2,length(trialOutcome));

figure()

for jj = 1:3
  
  %index by trial outcome
  idx_jj = ismember(spkCorr.condition, {['Accurate' trialOutcome{jj}],['Fast' trialOutcome{jj}]});  
  spkCorr_jj = spkCorr(idx_jj,:);
  
  %ANOVA -- Factors (groups) Condition (2) by Epoch (4)
  valsGroupsTbl = table();
  valsGroupsTbl.yVals = spkCorr_jj.(windowTest);
  valsGroupsTbl.condition = regexprep(spkCorr_jj.condition,'Correct|Error.*',''); %regular expression replace
  valsGroupsTbl.epoch = spkCorr_jj.alignedName;
  stats_Anova{jj} = satAnova(valsGroupsTbl);
  
  %prepare r_sc for plotting (index by epoch and task condition)
  idx_Acc = ismember(spkCorr_jj.condition, ['Accurate' trialOutcome{jj}]);
  idx_Fast = ismember(spkCorr_jj.condition, ['Fast' trialOutcome{jj}]);
  
  idx_Baseline = ismember(spkCorr_jj.alignedName, 'Baseline');
  idx_Visual = ismember(spkCorr_jj.alignedName, 'Visual');
  idx_PostSacc = ismember(spkCorr_jj.alignedName, 'PostSaccade');
  idx_PostRew = ismember(spkCorr_jj.alignedName, 'PostReward');
  
  %collect r_sc for plotting
  rsc_jj = spkCorr_jj.(windowTest);
  rsc_Acc{jj} = [rsc_jj(idx_Acc & idx_Baseline) , rsc_jj(idx_Acc & idx_Visual) , rsc_jj(idx_Acc & idx_PostSacc) , rsc_jj(idx_Acc & idx_PostRew)];
  rsc_Fast{jj} = [rsc_jj(idx_Fast & idx_Baseline) , rsc_jj(idx_Fast & idx_Visual) , rsc_jj(idx_Fast & idx_PostSacc) , rsc_jj(idx_Fast & idx_PostRew)];
  
  %compute mean and s.e.
  rsc_Acc_mu = nanmean(rsc_Acc{jj});
  rsc_Fast_mu = nanmean(rsc_Fast{jj});
  rsc_Acc_se = nanstd(rsc_Acc{jj}) / sqrt(size(rsc_Acc{jj},1));
  rsc_Fast_se = nanstd(rsc_Fast{jj}) / sqrt(size(rsc_Fast{jj},1));
  
  rsc_Baseline_Acc(:,jj) = [rsc_Acc_mu(1); rsc_Acc_se(1)];
  rsc_Baseline_Fast(:,jj) = [rsc_Fast_mu(1); rsc_Fast_se(1)];
  
  subplot(3,2,2*jj-1) %Accurate
  errorbar((1:4), rsc_Acc_mu, rsc_Acc_se, 'CapSize',0, 'Color','r')
  xticks(1:4); xticklabels([]); xlim([.8 4.2]); ylim(ylimPlot); ytickformat('%3.2f')
  
  subplot(3,2,2*jj) %Fast
  errorbar((1:4), rsc_Fast_mu, rsc_Fast_se, 'CapSize',0, 'Color',[0 .7 0])
  xticks(1:4); xticklabels([]); xlim([.8 4.2]); ylim(ylimPlot); yticks([])
  
  pause(0.25)
  
end % for : trialOutcome (jj)

title(areaPair)
ppretty([6.4,4])

%summary barplot
figure(); hold on
bar((1:2:5), rsc_Baseline_Acc(1,:), 'FaceColor','r', 'BarWidth',0.35)
errorbar((1:2:5), rsc_Baseline_Acc(1,:), rsc_Baseline_Acc(2,:), 'CapSize',0, 'Color','k')
bar((2:2:6), rsc_Baseline_Fast(1,:), 'FaceColor',[0 .7 0], 'BarWidth',0.35)
errorbar((2:2:6), rsc_Baseline_Fast(1,:), rsc_Baseline_Fast(2,:), 'CapSize',0, 'Color','k')
xticks([]); ytickformat('%3.2f')
ppretty([3,3])

clear idx_* rsc_* *jj

%% Save data
%save(outFn,'conditionByEpoch','pairAreasByEpoch');

%% Recode groups/factors for anova - PAIRAREA (3) by EPOCH (4) = (12*11)/2 = 66 comparisions
% valsGroupsTbl = table();
% valsGroupsTbl.yVals = spkCorrFilt.rhoRaw_200ms;
% valsGroupsTbl.pairAreas = spkCorrFilt.pairAreas;
% valsGroupsTbl.epoch = spkCorrFilt.alignedName;
% 
% [pairAreasByEpoch] = satAnova(valsGroupsTbl);





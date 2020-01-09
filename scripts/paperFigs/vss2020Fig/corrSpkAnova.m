%corrSpkAnova.m
% This script collects static rSC values for four trial epochs: Baseline,
% Visual Response, Post-Saccade, and Post-Reward. It computes an analysis
% of variance with factors Trial Epoch (4 levels) and Task Condition (2
% levels).

plot_X_epoch = false;
rhoColName = 'rhoRaw_150ms';
trialOutcome = {'Correct','ErrorChoice','ErrorTiming'};
lineStyle = {'-','--',':'};
ylimPlot = [0.04 0.18];

%load r_sc data
load('dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticNew.mat', 'SEF_FEF','SEF_SC','SEF_SEF');

%% Aggregate required fields into a single table
% pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
% Prune to columns of interest for each pairArea
useCols = {'pairAreas','XY_Dist','condition','alignedName',rhoColName};

%isolate area pairing of interest
% spkCorr = SEF_SEF(:,useCols);
spkCorr = [SEF_FEF(:,useCols); SEF_SC(:,useCols)];
areaPair = 'SEF-FEF & SEF-SC';

%take absolute value of correlation
spkCorr.(rhoColName) = abs(spkCorr.(rhoColName));

%initialize output for plot and ANOVA
rsc_Acc = cell(1,length(trialOutcome));
rsc_Fast = cell(1,length(trialOutcome));

rsc_Epoch_Acc = NaN(2,length(trialOutcome));
rsc_Epoch_Fast = NaN(2,length(trialOutcome));

if (plot_X_epoch)
  figure()
end

for jj = 1:3
  
  %index by trial outcome
  idx_jj = ismember(spkCorr.condition, {['Accurate' trialOutcome{jj}],['Fast' trialOutcome{jj}]});  
  spkCorr_jj = spkCorr(idx_jj,:);
  
  %prepare r_sc for plotting (index by epoch and task condition)
  idx_Acc = ismember(spkCorr_jj.condition, ['Accurate' trialOutcome{jj}]);
  idx_Fast = ismember(spkCorr_jj.condition, ['Fast' trialOutcome{jj}]);
  
  idx_Baseline = ismember(spkCorr_jj.alignedName, 'Baseline');
  idx_Visual = ismember(spkCorr_jj.alignedName, 'Visual');
  idx_PostSacc = ismember(spkCorr_jj.alignedName, 'PostSaccade');
  idx_PostRew = ismember(spkCorr_jj.alignedName, 'PostReward');
  
  %collect r_sc for plotting
  rsc_jj = spkCorr_jj.(rhoColName);
  rsc_Acc{jj} = [rsc_jj(idx_Acc & idx_Baseline) , rsc_jj(idx_Acc & idx_Visual) , rsc_jj(idx_Acc & idx_PostSacc) , rsc_jj(idx_Acc & idx_PostRew)];
  rsc_Fast{jj} = [rsc_jj(idx_Fast & idx_Baseline) , rsc_jj(idx_Fast & idx_Visual) , rsc_jj(idx_Fast & idx_PostSacc) , rsc_jj(idx_Fast & idx_PostRew)];
  
  %compute mean and s.e.
  rsc_Acc_mu = nanmean(rsc_Acc{jj});
  rsc_Fast_mu = nanmean(rsc_Fast{jj});
  rsc_Acc_se = nanstd(rsc_Acc{jj}) / sqrt(size(rsc_Acc{jj},1));
  rsc_Fast_se = nanstd(rsc_Fast{jj}) / sqrt(size(rsc_Fast{jj},1));
  
  %isolate rsc values for a particular epoch for barplot
  %1=Baseline, 2=VisResponse, 3=PostSacc, 4=PostRew
  epoch_Plot = 1; epoch_Print = 'Baseline';
  rsc_Epoch_Acc(:,jj) = [rsc_Acc_mu(epoch_Plot); rsc_Acc_se(epoch_Plot)];
  rsc_Epoch_Fast(:,jj) = [rsc_Fast_mu(epoch_Plot); rsc_Fast_se(epoch_Plot)];
  
  if (plot_X_epoch)
    subplot(1,2,1); hold on; title(areaPair) %Accurate
    errorbar((1:4), rsc_Acc_mu, rsc_Acc_se, 'CapSize',0, 'Color','r', 'LineStyle',lineStyle{jj})
    xticks(1:4); xticklabels([]); xlim([.8 4.2]); ylim(ylimPlot); ytickformat('%3.2f')
    subplot(1,2,2); hold on; title(areaPair) %Fast
    errorbar((1:4), rsc_Fast_mu, rsc_Fast_se, 'CapSize',0, 'Color',[0 .7 0], 'LineStyle',lineStyle{jj})
    xticks(1:4); xticklabels([]); xlim([.8 4.2]); ylim(ylimPlot); yticks([])
  end
  
  pause(0.25)
  
end % for : trialOutcome (jj)

if (plot_X_epoch)
  ppretty([6.4,1.4])
end

%summary barplot
figure(); hold on; title([areaPair, ' -- ', epoch_Print])
bar((1:2:5),      rsc_Epoch_Acc(1,:), 'FaceColor','r', 'BarWidth',0.35)
errorbar((1:2:5), rsc_Epoch_Acc(1,:), rsc_Epoch_Acc(2,:), 'CapSize',0, 'Color','k')
bar((2:2:6),      rsc_Epoch_Fast(1,:), 'FaceColor',[0 .7 0], 'BarWidth',0.35)
errorbar((2:2:6), rsc_Epoch_Fast(1,:), rsc_Epoch_Fast(2,:), 'CapSize',0, 'Color','k')
xticks([]); ytickformat('%3.2f')
ppretty([3,3])

%% Anova computation for different epochs etc...
%Two-way between-subjects ANOVA with factors Condition (Fast, Accurate) and
%Trial Outcome (Correct, Choice Error, Timing Error)
epochs = {'Baseline','Visual','PostSaccade','PostReward'};
rhoColName = 'rhoRaw_150ms';
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
statsAnova = struct();
for ep = 1:numel(epochs)
    epoch = epochs{ep};
    idx = ismember(spkCorr.alignedName,epoch);
    anovaTbl = table();
    anovaTbl.yVals = spkCorr{idx,rhoColName};
    anovaTbl.satCondition = spkCorr{idx,'satCondition'};
    anovaTbl.outcome = spkCorr{idx,'outcome'};
    statsAnova.(epoch) = satAnova(anovaTbl);
end

idx = ismember(spkCorr.condition, {['Accurate' trialOutcome{3}],['Fast' trialOutcome{3}]});
tempSpkCorrTbl = spkCorr(idx,:);
anovaTbl = table();
anovaTbl.yVals = spkCorr.(rhoColName);
anovaTbl.condition = regexprep(spkCorr_jj.condition,'Correct|Error.*',''); %regular expression replace
anovaTbl.epoch = spkCorr_jj.alignedName;
statsAnova{jj} = satAnova(anovaTbl);

clearvars -except SEF_* spkCorr stats_Anova rsc_Acc rsc_Fast

%% Save data
%save(outFn,'conditionByEpoch','pairAreasByEpoch');

%% Recode groups/factors for anova - PAIRAREA (3) by EPOCH (4) = (12*11)/2 = 66 comparisions
% valsGroupsTbl = table();
% valsGroupsTbl.yVals = spkCorrFilt.rhoRaw_200ms;
% valsGroupsTbl.pairAreas = spkCorrFilt.pairAreas;
% valsGroupsTbl.epoch = spkCorrFilt.alignedName;
% 
% [pairAreasByEpoch] = satAnova(valsGroupsTbl);




% % %summary barplot
% % epoch_Plot = 1; %1=Baseline, 2=VisResponse, 3=PostSacc, 4=PostRew
% % 
% % rsc_Correct = [rsc_Acc{1}(:,epoch_Plot), rsc_Fast{1}(:,epoch_Plot)];
% % num_Pair = size(rsc_Correct, 1);
% % rsc_ErrorChoice = [rsc_Acc{2}(:,epoch_Plot), rsc_Fast{2}(:,epoch_Plot)];
% % rsc_ErrorTiming = [rsc_Acc{3}(:,epoch_Plot), rsc_Fast{3}(:,epoch_Plot)];
% % 
% % figure(); hold on
% % bar([1 2], mean(rsc_Correct), 'FaceColor','b', 'BarWidth',0.35)
% % errorbar([1 2], mean(rsc_Correct), std(rsc_Correct)/sqrt(num_Pair), 'CapSize',0, 'Color','k')
% % bar([3 4], mean(rsc_ErrorChoice), 'FaceColor','b', 'BarWidth',0.35)
% % errorbar([3 4], mean(rsc_ErrorChoice), std(rsc_ErrorChoice)/sqrt(num_Pair), 'CapSize',0, 'Color','k')
% % bar([5 6], mean(rsc_ErrorTiming), 'FaceColor','b', 'BarWidth',0.35)
% % errorbar([5 6], mean(rsc_ErrorTiming), std(rsc_ErrorTiming)/sqrt(num_Pair), 'CapSize',0, 'Color','k')
% % xticks([]); ytickformat('%3.2f')
% % ppretty([3,3])


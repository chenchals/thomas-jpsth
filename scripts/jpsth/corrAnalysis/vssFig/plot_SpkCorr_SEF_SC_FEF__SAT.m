%plot_SpkCorr_SEF_SC_FEF__SAT.m

if (false)
  fileDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/';
  fileName = 'spkCorrAllPairsStaticNew.mat';
  load([fileDir fileName], 'SEF_FEF','SEF_SC','SEF_SEF')
end

PAIR = {spkCorr.SEF_SEF, spkCorr.SEF_FEF, spkCorr.SEF_SC};
pairName = {'SEF-SEF','SEF-FEF','SEF-SC'};
LINESTYLE = {'-', '--', ':'};
OFFSET = [-.03, 0, .03];

%combine SC and FEF
PAIR{2} = [PAIR{2}; PAIR{3}];   PAIR{3} = [];
pairName{2} = 'SEF-FEF/SC';     pairName{3} = [];

FIELD_RHO = 'rhoRaw_200ms';
FIELD_PVAL = 'pvalRaw_200ms';
% conds = {'Correct','ErrorChoice','ErrorTiming'};
conds = {'ErrorChoice'};

for c = 1:numel(conds)
  cond = conds{c};
  
  figure(); hold on
  
  doAbs = 1;
  drop0Dist = 1;
  
  for pp = 1:2

    %index by epoch
    idx_Baseline = ismember(PAIR{pp}.alignedName, 'Baseline');
    idx_Visual = ismember(PAIR{pp}.alignedName, 'Visual');
    idx_PostSacc = ismember(PAIR{pp}.alignedName, 'PostSaccade');
    idx_PostRew = ismember(PAIR{pp}.alignedName, 'PostReward');

    %index by condition and trial outcome
  %   idx_AccCorr = ismember(PAIR{pp}.condition, 'AccurateCorrect');
  %   idx_FastCorr = ismember(PAIR{pp}.condition, 'FastCorrect');
  %   idx_AccErrChc = ismember(PAIR{pp}.condition, 'AccurateErrorChoice');
  %   idx_FastErrChc = ismember(PAIR{pp}.condition, 'FastErrorChoice');
    idx_Acc = ismember(PAIR{pp}.condition, ['Accurate' cond]);
    idx_Fast = ismember(PAIR{pp}.condition, ['Fast' cond]);

    if drop0Dist
        idx_Acc = idx_Acc &  PAIR{pp}.XY_Dist ~= 0 ;
        idx_Fast = idx_Fast &  PAIR{pp}.XY_Dist ~= 0;     
    end

    %get spike count corr. coeff.
    rsc_Acc = [PAIR{pp}.(FIELD_RHO)(idx_Acc & idx_Baseline), PAIR{pp}.(FIELD_RHO)(idx_Acc & idx_Visual), ...
        PAIR{pp}.(FIELD_RHO)(idx_Acc & idx_PostSacc), PAIR{pp}.(FIELD_RHO)(idx_Acc & idx_PostRew)];
    rsc_Fast = [PAIR{pp}.(FIELD_RHO)(idx_Fast & idx_Baseline), PAIR{pp}.(FIELD_RHO)(idx_Fast & idx_Visual), ...
        PAIR{pp}.(FIELD_RHO)(idx_Fast & idx_PostSacc), PAIR{pp}.(FIELD_RHO)(idx_Fast & idx_PostRew)];

    if doAbs
        rsc_Acc = abs(rsc_Acc);
        rsc_Fast = abs(rsc_Fast);
    end


    %compute mean and se
    rsc_mu_Acc = nanmean(rsc_Acc);
    rsc_mu_Fast = nanmean(rsc_Fast);
    rsc_se_Acc = nanstd(rsc_Acc) / sqrt(size(rsc_Acc,1));
    rsc_se_Fast = nanstd(rsc_Fast) / sqrt(size(rsc_Fast,1));

    %plotting
    errorbar((1:4)+OFFSET(pp), rsc_mu_Acc, rsc_se_Acc, 'CapSize',0, 'Color','r', 'LineStyle',LINESTYLE{pp})
    errorbar((1:4)-OFFSET(pp), rsc_mu_Fast, rsc_se_Fast, 'CapSize',0, 'Color',[0 .7 0], 'LineStyle',LINESTYLE{pp})

    if (false)
    %stats -- ANOVA
    fprintf([pairName{pp}, '\n'])
    X = [rsc_Acc ; rsc_Fast];
    anova2_TR(X, size(rsc_Acc,1), 1, 3)
    end

  end % for : pair (pp)
  
  legend({'SEF-SEF','','SEF-FEF/SC',''}, 'Location','Northwest')
  xticks(1:4); xticklabels({'Baseline','Visual','PostSacc','PostRew'})
  title(cond)
  ppretty([4.8,2.4])
  ylim([0.0 0.25])
  
end % for : condition (c)

clearvars -except SEF_FEF SEF_SC SEF_SEF spkCorr




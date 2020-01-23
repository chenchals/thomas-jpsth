% Create different types of SDF heatmaps 
% Data that needs to be loaded are for all pairs from
% summary/spkCorrAllPairsStaticRhoPval.mat 
%% Load static RSC matrix for all pairs of SEF units
    fn = 'dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticRhoPval.mat';
        % Load data
    temp = load(fn,'-regexp','SEF*');
    allSpkCorr = table();
    fns = fieldnames(temp);
    for jj = 1:numel(fns)
        allSpkCorr = [allSpkCorr;temp.(fns{jj})]; %#ok<*AGROW>
    end
    clearvars temp fns jj
%% Plot all units by significance noFrOrSign filter
tic
unitsNoFilterTbl = categorizeUnitsByRscSignif(allSpkCorr);
useNormalized = 0; subDirName = 'noFilterMean';
[sdfsNoFilterMean,sdfsImgsNoFilterMean] = plotHeatmapsNoFilter(unitsNoFilterTbl,useNormalized,subDirName);
useNormalized = 1; subDirName = 'noFilterNorm';
[sdfsNoFilter,sdfsImgsNoFilter] = plotHeatmapsNoFilter(unitsNoFilterTbl,useNormalized,subDirName);
toc
%%

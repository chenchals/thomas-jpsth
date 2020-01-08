%% Jan 06, 2020
% 2.       What types of neurons contribute to significant r_sc at any point during the trial?
%   a.       Summary plot with SDFs of all such neurons in SEF.
%   b.       Summary plot for FEF.
%   c.       Summary plot for SC.
%   d.       Thomas & Chenchal to meet in 069 to discuss individual SDFs after viewing summary plot.
% 3.       What types of neurons do not contribute to significant r_sc at any point during the trial?
%   a.       Summary plot with SDFs of all such neurons in SEF.
%   b.       Summary plot for FEF.
%   c.       Summary plot forSC.
%   d.       Thomas & Chenchal to meet in 069.
%  
% Yes, and let?s also implement the color density SDF raster plot so we can see variation across all neurons in given groups
%  
%%
unitsTbl = categorizeUnitsByRscSignif();
% Filter criteria used to categorize units by Rsc
filterEpochs = unique(unitsTbl.filter_Epoch);
filterOutcomes = unique(unitsTbl.filter_Outcome);
filterPval = unique(unitsTbl.filter_Pval);

satSdfDir = 'dataProcessed/dataset/satSdfs';

%% Use results from:
epoch = 'PostSaccade';
outcome = 'Correct';
pval = 0.05;
currUnitsByArea = getUnitNums(unitsTbl,epoch,outcome,pval);
unitAreas = {'SEF','FEF','SC'}; % sameArea_SEF_X, sameArea_SEF_Y
sdfsTbl = table();
isSignifs = [1 0];
for s = 1:numel(isSignifs)
    signif = isSignifs(s);
    idx = currUnitsByArea.filter_IsRscSignificant == signif;
    for a = 1:numel(unitAreas)
        unitNums = currUnitsByArea.SEF{idx};
        parfor un = 1:numel(unitNums)
            % compute average SDF for unit by condition and epoch
            unitNum = unitNums(un);
            temp = getNormalizedSdf(satSdfDir,unitNum);
            temp.isRscSignificant = repmat(signif,size(temp,1),1);
            sdfsTbl = [sdfsTbl;temp];
        end
    end
end


%% Other functions
function [sigNonSigUnits] = getUnitNums(unitsTbl,epoch,outcome,pval)
    idx = ismember(unitsTbl.filter_Epoch,epoch) & ismember(unitsTbl.filter_Outcome,outcome) & unitsTbl.filter_Pval == pval;
    sigNonSigUnits = unitsTbl(idx,{'filter_IsRscSignificant','SEF','FEF','SC'}) ;
end

function [satSdfTbl] = getNormalizedSdf(satSdfDir,unitNum)
    satSdfTbl = table();
    baselineWin = [-500 -100];
    fn = fullfile(satSdfDir,num2str(unitNum,'Unit_%03d.mat'));
    sdfs = load(fn,'sdfs');
    sdfs = sdfs.sdfs;
    % timestamps for SDFs
    tsV = sdfs.Visual_timeMs{1};
    % get Baseline FR - use visual epoch from all outcomes
    meanV = mean(cat(1,sdfs.Visual_sdfByTrial{:}),1);
    frBl = mean(meanV(find(tsV==baselineWin(1)):find(tsV==baselineWin(2))));
    % get max FR - use all epochs from all outcomes
    allFr = [meanV,...
             mean(cat(1,sdfs.PostSaccade_sdfByTrial{:}),1),...
             mean(cat(1,sdfs.PostReward_sdfByTrial{:}),1)];
    frMax = max(abs(allFr));
    % Create output table
    satSdfTbl.unitNum = sdfs.unitNum;
    satSdfTbl.condition = sdfs.condition;
    % do normalization for trail sdfs for all epochs for all outcomes
    % (trial-sdf - frBL)/(frMax-frBl)
    % Visual
    satSdfTbl.VisualAlignedEvent = sdfs.Visual_alignedEvent;
    satSdfTbl.VisualTs = sdfs.Visual_timeMs;
    t = arrayfun(@(x) (sdfs.Visual_sdfByTrial{x}-frBl)./(frMax-frBl),(1:size(sdfs,1))','UniformOutput',false);
    satSdfTbl.VisualSdf = cellfun(@(x) mean(x,1),t,'UniformOutput',false);
    % PostSaccade
    satSdfTbl.PostSaccadeAlignedEvent = sdfs.PostSaccade_alignedEvent;
    satSdfTbl.PostSaccadeTs = sdfs.PostSaccade_timeMs;
    t = arrayfun(@(x) (sdfs.PostSaccade_sdfByTrial{x}-frBl)./(frMax-frBl),(1:size(sdfs,1))','UniformOutput',false);
    satSdfTbl.PostSaccadeSdf = cellfun(@(x) mean(x,1),t,'UniformOutput',false);
    % PostReward
    satSdfTbl.PostRewardAlignedEvent = sdfs.PostReward_alignedEvent;
    satSdfTbl.PostRewardTs = sdfs.PostReward_timeMs;
    t = arrayfun(@(x) (sdfs.PostReward_sdfByTrial{x}-frBl)./(frMax-frBl),(1:size(sdfs,1))','UniformOutput',false);
    satSdfTbl.PostRewardSdf = cellfun(@(x) mean(x,1),t,'UniformOutput',false);

end

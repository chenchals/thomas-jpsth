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
parpoolSize = 0;
if ~isempty(gcp('nocreate'))
    parpoolSize = Inf;
end
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
        unitArea = unitAreas{a};
        unitNums = currUnitsByArea.(unitArea){idx};
        fprintf('Doing area %s...',unitArea);
        if isempty(unitNums)
            fprintf(' no units! Done!\n');
            continue;
        end
        parfor (un = 1:numel(unitNums), parpoolSize)
            % compute average SDF for unit by condition and epoch
            unitNum = unitNums(un);
            temp = getNormalizedSdf(satSdfDir,unitNum);
            temp.isRscSignificant = repmat(signif,size(temp,1),1);
            temp.unitArea = repmat({unitArea},size(temp,1),1);
            sdfsTbl = [sdfsTbl;temp];            
        end
        fprintf('Done!\n');
    end
end

%% Aggregate SDFs of all units by area by significance
sdfsImageTbl = table();
for s = 1:numel(isSignifs)
    signif = isSignifs(s);
    for a = 1:numel(unitAreas)
        unitArea = unitAreas{a};
        idx = sdfsTbl.isRscSignificant == signif ...
            & ismember(sdfsTbl.unitArea,unitArea); 
        conds = unique(sdfsTbl.condition,'stable');
        fprintf('Doing area %s...',unitArea);
        for c = 1:numel(conds)
            % aggregate all SDFs by condition for all epochs
            temp = table();
            condition = conds{c};
            cIdx = find(idx & ismember(sdfsTbl.condition,condition));
            if isempty(cIdx)
                fprintf(' no units! Done!\n');
                continue;
            end
            %cat(1,sdfs.Visual_sdfByTrial{:})
            temp.condition =  {condition};
            temp.isRscSignificant = signif;
            temp.unitArea = {unitArea};
            % Visual (includes Baseline)
            temp.VisualAlignedEvent = sdfsTbl.VisualAlignedEvent(cIdx(1));
            temp.VisualTs = sdfsTbl.VisualTs(cIdx(1));
            temp.VisualSdfs = {cat(1,sdfsTbl.VisualSdf{cIdx})};
            % PostSaccade 
            temp.PostSaccadeAlignedEvent = sdfsTbl.PostSaccadeAlignedEvent(cIdx(1));
            temp.PostSaccadeTs = sdfsTbl.PostSaccadeTs(cIdx(1));
            temp.PostSaccadeSdfs = {cat(1,sdfsTbl.PostSaccadeSdf{cIdx})};
            % PostReward 
            temp.PostRewardAlignedEvent = sdfsTbl.PostRewardAlignedEvent(cIdx(1));
            temp.PostRewardTs = sdfsTbl.PostRewardTs(cIdx(1));
            temp.PostRewardSdfs = {cat(1,sdfsTbl.PostRewardSdf{cIdx})};
            % Add row for condition
            sdfsImageTbl = [sdfsImageTbl;temp];            
        end
        fprintf('Done!\n');
    end
end

%% explore

%% Other functions
function [sigNonSigUnits] = getUnitNums(unitsTbl,epoch,outcome,pval)
    idx = ismember(unitsTbl.filter_Epoch,epoch) & ismember(unitsTbl.filter_Outcome,outcome) & unitsTbl.filter_Pval == pval;
    sigNonSigUnits = unitsTbl(idx,{'filter_IsRscSignificant','SEF','FEF','SC'}) ;
end

function [satSdfTbl] = getNormalizedSdf(satSdfDir,unitNum)
    satSdfTbl = table();
    nTrialsForCondition = 50;
    baselineWin = [-600 0];
    fn = fullfile(satSdfDir,num2str(unitNum,'Unit_%03d.mat'));
    sdfs = load(fn,'sdfs');
    sdfs = sdfs.sdfs;
    % sdfs.outcome
    sdfs.outcome = regexprep(sdfs.condition,{'Fast','Accurate'},{'',''});
    % compute meanSdfs for all conditions by epoch
    sdfs.meanV = cellfun(@mean,sdfs.Visual_sdfByTrial,'UniformOutput',false);
    sdfs.meanPs = cellfun(@mean,sdfs.PostSaccade_sdfByTrial,'UniformOutput',false);
    sdfs.meanPr = cellfun(@mean,sdfs.PostReward_sdfByTrial,'UniformOutput',false);
    % 
    
    
    
    
    
    % timestamps for SDFs
    tsV = sdfs.Visual_timeMs{1};

    % 
    
    % get Baseline FR - use visual epoch from all outcomes
    meanV = mean(cat(1,sdfs.Visual_sdfByTrial{:}),1);
    frBl = mean(meanV(find(tsV==baselineWin(1)):find(tsV==baselineWin(2))));
    % get max FR - use all epochs from all outcomes
    meanPs = mean(cat(1,sdfs.PostSaccade_sdfByTrial{:}),1);
    meanPr = mean(cat(1,sdfs.PostReward_sdfByTrial{:}),1);
    meanAll = [meanV,meanPs,meanPr];
    frMaxAll = max(abs(meanAll));
    frMaxV = max(max(cat(1,sdfs.Visual_sdfByTrial{:})-frBl)); % there are no negative FRs
    frMaxPs = max(max(cat(1,sdfs.PostSaccade_sdfByTrial{:})-frBl)); % there are no negative FRs
    frMaxPr = max(max(cat(1,sdfs.PostReward_sdfByTrial{:})-frBl)); % there are no negative FRs
    % Create output table
    satSdfTbl.unitNum = sdfs.unitNum;
    satSdfTbl.condition = sdfs.condition;
    % do normalization for trail sdfs for all epochs for all outcomes
    % (trial-sdf - frBL)/(frMax-frBl)
    % Visual
    satSdfTbl.VisualAlignedEvent = sdfs.Visual_alignedEvent;
    satSdfTbl.VisualTs = sdfs.Visual_timeMs;
    t = arrayfun(@(x) (sdfs.Visual_sdfByTrial{x}-frBl)./(frMaxV),(1:size(sdfs,1))','UniformOutput',false);
    satSdfTbl.VisualSdf = cellfun(@(x) mean(x,1),t,'UniformOutput',false);
    % PostSaccade
    satSdfTbl.PostSaccadeAlignedEvent = sdfs.PostSaccade_alignedEvent;
    satSdfTbl.PostSaccadeTs = sdfs.PostSaccade_timeMs;
    t = arrayfun(@(x) (sdfs.PostSaccade_sdfByTrial{x}-frBl)./(frMaxPs),(1:size(sdfs,1))','UniformOutput',false);
    satSdfTbl.PostSaccadeSdf = cellfun(@(x) mean(x,1),t,'UniformOutput',false);
    % PostReward
    satSdfTbl.PostRewardAlignedEvent = sdfs.PostReward_alignedEvent;
    satSdfTbl.PostRewardTs = sdfs.PostReward_timeMs;
    t = arrayfun(@(x) (sdfs.PostReward_sdfByTrial{x}-frBl)./(frMaxPr),(1:size(sdfs,1))','UniformOutput',false);
    satSdfTbl.PostRewardSdf = cellfun(@(x) mean(x,1),t,'UniformOutput',false);

end

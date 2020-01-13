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
satSdfsTbl = table();
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
        %for un = 1:numel(unitNums)
            % compute average SDF for unit by condition and epoch
            unitNum = unitNums(un);
            temp = getNormalizedSdf(satSdfDir,unitNum);
            temp.isRscSignificant = repmat(signif,size(temp,1),1);
            temp.unitArea = repmat({unitArea},size(temp,1),1);
            satSdfsTbl = [satSdfsTbl;temp];            
        end
        fprintf('Done!\n');
    end
end

%% Aggregate SDFs of all units by area by significance
satSdfsImageTbl = table();
for s = 1:numel(isSignifs)
    signif = isSignifs(s);
    for a = 1:numel(unitAreas)
        unitArea = unitAreas{a};
        idx = satSdfsTbl.isRscSignificant == signif ...
            & ismember(satSdfsTbl.unitArea,unitArea); 
        conds = unique(satSdfsTbl.condition,'stable');
        fprintf('Doing area %s...',unitArea);
        for c = 1:numel(conds)
            % aggregate all SDFs by condition for all epochs
            temp = table();
            condition = conds{c};
            cIdx = find(idx & ismember(satSdfsTbl.condition,condition));
            if isempty(cIdx)
                fprintf(' no units! Done!\n');
                continue;
            end
            %cat(1,sdfs.Visual_sdfByTrial{:})
            temp.condition =  {condition};
            temp.isRscSignificant = signif;
            temp.unitArea = {unitArea};
            % Visual (includes Baseline)
            temp.VisualAlignedEvent = satSdfsTbl.VisualAlignedEvent(cIdx(1));
            temp.VisualTs = satSdfsTbl.VisualTs(cIdx(1));
            temp.VisualSdfs = {cat(1,satSdfsTbl.VisualSdf{cIdx})};
            % PostSaccade 
            temp.PostSaccadeAlignedEvent = satSdfsTbl.PostSaccadeAlignedEvent(cIdx(1));
            temp.PostSaccadeTs = satSdfsTbl.PostSaccadeTs(cIdx(1));
            temp.PostSaccadeSdfs = {cat(1,satSdfsTbl.PostSaccadeSdf{cIdx})};
            % PostReward 
            temp.PostRewardAlignedEvent = satSdfsTbl.PostRewardAlignedEvent(cIdx(1));
            temp.PostRewardTs = satSdfsTbl.PostRewardTs(cIdx(1));
            temp.PostRewardSdfs = {cat(1,satSdfsTbl.PostRewardSdf{cIdx})};
            % Add row for condition
            satSdfsImageTbl = [satSdfsImageTbl;temp];            
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

function [satSdfsTbl] = getNormalizedSdf(satSdfDir,unitNum)
    baselineWin = [-600 0];
    fn = fullfile(satSdfDir,num2str(unitNum,'Unit_%03d.mat'));
    sdfs = load(fn,'sdfs');
    sdfs = sdfs.sdfs;
    % sdfs.satCondition
    sdfs.satCondition = regexprep(sdfs.condition,{'Correct','Error.*'},{'',''});
    % sdfs.outcome
    sdfs.outcome = regexprep(sdfs.condition,{'Fast','Accurate'},{'',''});
    % timestamps for SDFs
    tsV = sdfs.Visual_timeMs{1};
    % get Baseline FR - use visual epoch from all conditions (6 conditions)
    meanV = mean(cat(1,sdfs.Visual_sdfByTrial{:}),1);
    frBl = mean(meanV(find(tsV==baselineWin(1)):find(tsV==baselineWin(2))));
    % get max FR - use all epochs from all outcomes
    meanPs = mean(cat(1,sdfs.PostSaccade_sdfByTrial{:}),1);
    meanPr = mean(cat(1,sdfs.PostReward_sdfByTrial{:}),1);
    meanAll = [meanV,meanPs,meanPr];
    frMaxAll = max(abs(meanAll));
    % do normalization for trail sdfs for all outcomes by SatCondition
    % For Example Fast Visual =
    % all-FastVisual trials=[FastErrorChoice, FastErrorTiming, FastErrorChoice]
    % (allFastVisual-trial-sdfs - frBL)/(frMax-frBl)
    satConds = unique(sdfs.satCondition);
    epochs = {'Visual','PostSaccade','PostReward'};
    satSdfsTbl = table();    
    for sc = 1:numel(satConds)
        tempSatSdfTbl = table();   
        satCondition = satConds{sc};
        tempSatSdfTbl.unitNum = unitNum;
        tempSatSdfTbl.condition = {satCondition};
        tempSatSdfTbl.frBaseline = frBl;
        tempSatSdfTbl.frMaxAllConditionsEpochs = frMaxAll;       
        idxSatCond = ismember(sdfs.satCondition,satCondition);
        for ep = 1:numel(epochs)
            epoch = epochs{ep};
            meanSatSdf = mean(cat(1,sdfs.([epoch '_sdfByTrial']){idxSatCond}),1);
            normSatSdf = (meanSatSdf-frBl)./(frMaxAll-frBl);
            tempSatSdfTbl.([epoch 'AlignedEvent']) = sdfs.([epoch '_alignedEvent']){1};
            tempSatSdfTbl.([epoch 'TimeMs']) = sdfs.([epoch '_timeMs']){1};
            tempSatSdfTbl.([epoch 'SatSdfMean']) = meanSatSdf;
            tempSatSdfTbl.([epoch 'SatSdfNormalized']) = normSatSdf;
        end  
        satSdfsTbl = [satSdfsTbl;tempSatSdfTbl]; %#ok<*AGROW>
    end
end

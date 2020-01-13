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
filter.epoch = 'PostSaccade';
filter.outcome = 'Correct';
filter.pval = 0.05;
currUnitsByArea = getUnitNums(unitsTbl,filter.epoch,filter.outcome,filter.pval);
unitAreas = {'SEF','FEF','SC'}; % sameArea_SEF_X, sameArea_SEF_Y
satSdfsTbl = table();
isSignifs = [1 0];
for sig = 1:numel(isSignifs)
    signif = isSignifs(sig);
    idx = currUnitsByArea.filter_IsRscSignificant == signif;
    for ar = 1:numel(unitAreas)
        unitArea = unitAreas{ar};
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
            tempSatSdfTbl = getNormalizedSdf(satSdfDir,unitNum); 
            tempSatSdfTbl.isRscSignificant = repmat(signif,size(tempSatSdfTbl,1),1);
            tempSatSdfTbl.unitArea = repmat({unitArea},size(tempSatSdfTbl,1),1);
            satSdfsTbl = [satSdfsTbl;tempSatSdfTbl];            
        end
        fprintf('Done!\n');
    end
end

%% Aggregate SDFs of all units by (1)significance of r_sc, (2)area, (3)SAT condition (across all outcomes), and (4)epochs
satSdfsImageTbl = table();
satConds = {'Fast','Accurate'};
epochs = {'Visual','PostSaccade','PostReward'};
for sig = 1:numel(isSignifs)
    signif = isSignifs(sig);
    idxSig = satSdfsTbl.isRscSignificant == signif;
    for ar = 1:numel(unitAreas)
        unitArea = unitAreas{ar};
        fprintf('Aggregating SDFs for heatmap area %s...',unitArea);
        idxAr = ismember(satSdfsTbl.unitArea,unitArea);
        for sc = 1:numel(satConds)
            satCondition = satConds{sc};
            idxSat = ismember(satSdfsTbl.satCondition,satCondition);
            idx = find(idxSig & idxAr & idxSat);
            if isempty(idx)
                fprintf('[0 units...]');
                continue;
            end
            fprintf('[%d units...]',numel(idx));
            for ep = 1:numel(epochs)
                epoch = epochs{ep};
                tempTbl = table();
                % Filter used for finding significant pairs
                tempTbl.filterOutcome = {filter.outcome};
                tempTbl.filterEpoch = {filter.epoch};
                tempTbl.filterPval = {filter.pval};
                % 
                tempTbl.unitArea = {unitArea};
                tempTbl.unitNums = {satSdfsTbl.unitNum(idx)};
                tempTbl.numUnits = numel(tempTbl.unitNums{1});
                tempTbl.isRscSignificant = signif;
                tempTbl.satCondition = {satCondition};                
                tempTbl.frBaseline =  {satSdfsTbl.frBaseline(idx)};
                tempTbl.frMaxAllConditionsEpochs = {satSdfsTbl.frMaxAllConditionsEpochs(idx)};       
                tempTbl.([epoch 'AlignedEvent']) = satSdfsTbl.([epoch 'AlignedEvent'])(idx(1));
                tempTbl.([epoch 'TimeMs']) = satSdfsTbl.([epoch 'TimeMs'])(idx(1));
                tempTbl.([epoch 'SatSdfMean']) = {cat(1,satSdfsTbl.([epoch 'SatSdfMean']){idx})};
                tempTbl.([epoch 'SatSdfNormalized']) = {cat(1,satSdfsTbl.([epoch 'SatSdfNormalized']){idx})};
           
            end
            % Add row for condition
            satSdfsImageTbl = [satSdfsImageTbl;tempTbl];            
        end
        fprintf('Done!\n');
    end
end
%% Show satSdfsImageTbl




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
        tempSatSdfTbl.satCondition = {satCondition};
        tempSatSdfTbl.frBaseline = frBl;
        tempSatSdfTbl.frMaxAllConditionsEpochs = frMaxAll;       
        idxSatCond = ismember(sdfs.satCondition,satCondition);
        for ep = 1:numel(epochs)
            epoch = epochs{ep};
            meanSatSdf = mean(cat(1,sdfs.([epoch '_sdfByTrial']){idxSatCond}),1);
            normSatSdf = (meanSatSdf-frBl)./(frMaxAll-frBl);
            tempSatSdfTbl.([epoch 'AlignedEvent']) = sdfs.([epoch '_alignedEvent'])(1);
            tempSatSdfTbl.([epoch 'TimeMs']) = {sdfs.([epoch '_timeMs']){1}};
            tempSatSdfTbl.([epoch 'SatSdfMean']) = {meanSatSdf};
            tempSatSdfTbl.([epoch 'SatSdfNormalized']) = {normSatSdf};
        end  
        satSdfsTbl = [satSdfsTbl;tempSatSdfTbl]; %#ok<*AGROW>
    end
end

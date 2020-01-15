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
isSignifs = [1 0];
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
            tempTbl = table();
            for ep = 1:numel(epochs)
                epoch = epochs{ep};

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
[H_plots,H_Figure] = getPlotHandles();
isSignifs = [1 0];
satConds = {'Fast','Accurate'};
epochs = {'Visual','PostSaccade','PostReward'};
unitAreas = {'SEF', 'FEF','SC'};
plotNo = 0;
% use red to yellow (dark -> white) (0 -> 1)
% Each row(6 plots): [(signif),V, PS, PR], [(not-signif),V, PS, PR]
% Each col(6 plots): [Fast- SEF,FEF,SC], [Accurate- SEF, FEF, SC]
pltSatCond = 'Fast';
pltArea = 'SEF';
% Plots are plotted column-wise
for sig = 1:numel(isSignifs)
    signif = isSignifs(sig);
    for ep = 1:numel(epochs)
        epoch = epochs{ep};
        for sc = 1:numel(satConds)
            satCondition = satConds{sc};
            for ar = 1:numel(unitAreas)
                unitArea = unitAreas{ar};
                plotNo = plotNo + 1;
                idx = find(ismember(satSdfsImageTbl.satCondition,satCondition)...
                    & satSdfsImageTbl.isRscSignificant == signif ...
                    & ismember(satSdfsImageTbl.unitArea,unitArea));
                if isempty(idx)
                    % not units for the area...
                    continue;
                end
                tempSdf = satSdfsImageTbl(idx,:);
                % get the SDFs and related information for plot
                sdfImg = tempSdf.([epoch 'SatSdfNormalized']){1};
                timeMs = tempSdf.([epoch 'TimeMs']){1};
                yLabel = pltArea;
                showImage(H_plots(plotNo),sdfImg,timeMs,epoch,unitArea);
            end % for each area SEF,FEF,SC...
        end % for each SAT condition Fast, Accurate
    end % for each epoch Visual, PostSaccade, PostReward
end % for each significance level 1=isSignificant, 0=isNotSignificant
%% explore

%% Other functions
function [] = showImage(H_axis,sdfImg,timeMs,epoch,pltArea)
    axes(H_axis);
    yMax = size(sdfImg,1);
    % show image
    hImg = imagesc(sdfImg);
    % row 1 = top of image, ie the SDF of sdfImg(1,:)
    set(gca,'YDir','reverse');
    xTickLabel = unique([0:-200:min(timeMs) 0:200:max(timeMs)]);
    if ~xTickLabel(1) % first label is zero
        xTickLabel = xTickLabel(1:end-1);
    else
        xTickLabel = xTickLabel(2:end-1);
    end
    % location of ticks on the image to correspond to tick labels
    % xTickLoc = linspace(1, size(sdfImg, 2), numel(xTickLabel));
    xTickLoc = arrayfun(@(x) find(timeMs==x),xTickLabel);
    % Tick location for 0 ms (align Event time)
    x0Loc = xTickLoc(xTickLabel==0);
    % add zero line
    line([x0Loc x0Loc],[0 yMax+1]);
    
    set(gca, 'XTick', xTickLoc)
    set(gca,'XTickLabel', xTickLabel); 
    set(gca,'XMinorGrid','on')    

end

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

% Get figure template
function [H_plots,H_Figure] = getPlotHandles()
    % For Heatmap of aggregated SDFs for each area
    % 6 rows {3-Areas for Fast, 3-areas for Accurate}
    % 6 cols {3-epochs for significant, 3-epochs for non-significant}
    H_Figure = newFigure();
    nRows = 6;
    startOffset = 0.06;
    pltH = 0.13;
    pltWFor3Cols = 0.4;
    gutterx = 0.005;
    guttery = 0.015;
    % conditions SDFs have different lengths of times (x axis).
    % make widths proportional such that time unit ticks are of equal length
    % visual:[-600 400] , postSaccade:[-200 600], postReward:[-100 400]
    pltWProp = [1001 801 501];
    % partition total plot width to the proportion above
    pltWs = pltWProp.*(pltWFor3Cols/sum(pltWProp));
    offsetsX = [0 pltWs(1)+gutterx sum(pltWs(1:2))+gutterx*2] + startOffset;
    offsetsY = 0.95-pltH:-(pltH+guttery):guttery; % for 6 row-starts
    pltCount = 0;
    for col = 1:numel(offsetsY)
        if col < 4
            pltW = pltWs(col);
            pos(1) = offsetsX(col);
        else
            pltW = pltWs(col-3);
            pos(1) = (pltWFor3Cols + gutterx*16) + offsetsX(col-3);
        end
        pos(3:4) = [pltW pltH];
        for ro = 1:nRows
            if ro < 4
                pos(2) = offsetsY(ro);
            else
                pos(2) = offsetsY(ro) - guttery*2;
            end
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount));
        end
    end
end
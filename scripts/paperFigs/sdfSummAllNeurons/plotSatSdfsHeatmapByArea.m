function [satSdfsTbl, satSdfsImageTbl] = plotSatSdfsHeatmapByArea(unitsTbl,useOutcome,useEpoch,usePval,useNormalized,subDirName)

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
outputPdfDir = fullfile('dataProcessed/analysis/spkCorr/summary/sdfHeatmaps',subDirName);
if ~exist(outputPdfDir,'dir')
    mkdir(outputPdfDir);
end
satSdfDir = 'dataProcessed/dataset/satSdfs';

%% Use results from:
parpoolSize = 0;
if ~isempty(gcp('nocreate'))
    parpoolSize = Inf;
end
useFilter.epoch = useEpoch; %'PostSaccade';
useFilter.outcome = useOutcome; %'Correct';
useFilter.pval = usePval; %0.05;

fprintf('Doing plotSatSdfsHeatmapByArea for outcome = %s, epoch = %s, and pval <= %0.02f\n',...
        useFilter.outcome,useFilter.epoch,useFilter.pval)

currUnitsByArea = getUnitNums(unitsTbl,useFilter.epoch,useFilter.outcome,useFilter.pval);
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
        unitsSortOrder = [];
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
                tempTbl.filterOutcome = {useFilter.outcome};
                tempTbl.filterEpoch = {useFilter.epoch};
                tempTbl.filterPval = {useFilter.pval};
                %
                tempTbl.unitArea = {unitArea};
                tempTbl.unitNums = {satSdfsTbl.unitNum(idx)};
                tempTbl.numUnits = numel(tempTbl.unitNums{1});
                tempTbl.isRscSignificant = signif;
                tempTbl.satCondition = {satCondition};
                tempTbl.frBaseline =  {satSdfsTbl.frBaseline(idx)};
                tempTbl.frBaselineStd =  {satSdfsTbl.frBaselineStd(idx)};
                tempTbl.frMaxAllConditionsEpochs = {satSdfsTbl.frMaxAllConditionsEpochs(idx)};
                tempTbl.([epoch 'AlignedEvent']) = satSdfsTbl.([epoch 'AlignedEvent'])(idx(1));
                tempTbl.([epoch 'TimeMs']) = satSdfsTbl.([epoch 'TimeMs'])(idx(1));
                tempTbl.([epoch 'SatSdfMean']) = {cat(1,satSdfsTbl.([epoch 'SatSdfMean']){idx})};
                tempTbl.([epoch 'SatSdfNormalized']) = {cat(1,satSdfsTbl.([epoch 'SatSdfNormalized']){idx})}; 
                
                % add a UnitSortOrderByVisual for Visual epoch
                if strcmp(epoch,'Visual')
                    if sc == 1 % use only one SAT condition for sorting
                        visStartWin = find(tempTbl.VisualTimeMs{1}==0);
                        visEndWin = visStartWin + 300;
                        [~, unitsSortOrder] = sortSdfsMat(tempTbl.VisualSatSdfNormalized{1},visStartWin,visEndWin);
                    end
                    tempTbl.VisualSortOrderForUnits = {unitsSortOrder};
                end
            end
            % Add row for condition
            satSdfsImageTbl = [satSdfsImageTbl;tempTbl];
        end
        fprintf('Done!\n');
    end
end
%% For normalized plots set CLimits
cLims = [-1 1]; % so as to account for supression
if ~useNormalized
    % find grand men and max for all the sdfs across all areas,
    % epochs,signifs, sat
    allSdfMats = [satSdfsImageTbl.VisualSatSdfMean;...
                  satSdfsImageTbl.PostSaccadeSatSdfMean;...
                  satSdfsImageTbl.PostRewardSatSdfMean];
    cLims(1) = min(cellfun(@(x) min(min(x)),allSdfMats)); 
    cLims(2) = max(cellfun(@(x) max(max(x)),allSdfMats)); 
end

%% Show satSdfsImageTbl
% use satSdfsImageTbl
[H_plots,H_Figure] = getPlotHandles();
fClr = [0 0.7 0];
aClr = [1 0 0];
isSignifs = [1 0];
satConds = {'Fast','Accurate'};
epochs = {'Visual','PostSaccade','PostReward'};
alignedEvents = {'CueOn','SaccadePrimary','Reward'};
unitAreas = {'SEF', 'FEF','SC'};
plotNo = 0;
% use red to yellow (dark -> white) (0 -> 1)
% Each row(6 plots): [(signif),V, PS, PR], [(not-signif),V, PS, PR]
% Each col(6 plots): [Fast- SEF,FEF,SC], [Accurate- SEF, FEF, SC]
% Plots are plotted column-wise
for sig = 1:numel(isSignifs)
    signif = isSignifs(sig);
    annotateSignif = 1;
    for ep = 1:numel(epochs)
        epoch = epochs{ep};
        alignedEvent = alignedEvents{ep};
        for sc = 1:numel(satConds)
            satCondition = satConds{sc};
            annotateSat = 1;
            for ar = 1:numel(unitAreas)
                unitArea = unitAreas{ar};
                plotNo = plotNo + 1;
                idx = find(ismember(satSdfsImageTbl.satCondition,satCondition)...
                    & satSdfsImageTbl.isRscSignificant == signif ...
                    & ismember(satSdfsImageTbl.unitArea,unitArea));
                if isempty(idx)
                    % not units for the area...
                    %set(H_plots(plotNo),'Visible','off');
                    set(H_plots(plotNo),'XColor',[1 1 1],'YColor',[1 1 1],'YTickLabel',{});
                    h = get(H_plots(plotNo),'YLabel');
                    set(h,'String',{sprintf('%s ',unitArea),sprintf('[%d] ',0)});
                    set(h,'Rotation',0,'FontWeight','bold','HorizontalAlignment','right');
                    set(h,'Color','k');
                    continue;
                end
                % get the SDFs and related information for plot
                tempSdf = satSdfsImageTbl(idx,:);
                if useNormalized
                    sdfImg = tempSdf.([epoch 'SatSdfNormalized']){1};
                else
                    sdfImg = tempSdf.([epoch 'SatSdfMean']){1};
                end
                timeMs = tempSdf.([epoch 'TimeMs']){1};
                sortOrder = tempSdf.VisualSortOrderForUnits{1};
                sdfImg = sdfImg(sortOrder,:);
                showImage(H_plots(plotNo),sdfImg,timeMs,epoch,alignedEvent,unitArea,cLims);                
                % set xticklabel of this plot if exists
                if ar>1
                    set(H_plots(plotNo-1),'XTickLabel',{});
                    delete(get(H_plots(plotNo-1),'XLabel'));
                    %delete title from current plot
                    delete(get(H_plots(plotNo),'Title'));
                end
                % Annotate SAT condition
                if annotateSat
                    str = upper(satCondition);
                    clr = fClr;
                    pos = [0.01 0.91 0.3 0.04];
                    if sc == 2
                        pos = [0.01 0.45 0.3 0.04];
                        clr = aClr;
                    end
                    if sig == 2
                        pos(1) = 0.501;
                    end
                    annotation('textbox','String',str,'FontSize',18,'FontWeight','bold',...
                        'Color',clr,'Position',pos,'LineStyle','none');
                    annotateSat = 0;
                end
                % Annotate Significant vs non-significant for the top row
                if annotateSignif
                    str = upper('Significant Units');
                    pos = [0.15 0.925 0.3 0.04];
                    if signif == 0
                        str = upper('Non significant Units');
                        pos = [0.65 0.925 0.3 0.04];
                    end
                    annotation('textbox','String',str,'FontSize',14,'FontWeight','bold',...
                        'Position',pos,'LineStyle','none','HorizontalAlignment','center');
                    annotateSignif = 0;
                end
                
            end % for each area SEF,FEF,SC...
        end % for each SAT condition Fast, Accurate
    end % for each epoch Visual, PostSaccade, PostReward
end % for each significance level 1=isSignificant, 0=isNotSignificant
% Cleanup labels and ticks for the 6 by 3 and another 6 by 3 set of plots
pltNos = reshape(1:36,6,6);
% No YLabels for all except for 1 and 4 columns
pltIdx = pltNos(:,[2,3,5,6]);
arrayfun(@(x) delete(get(H_plots(x),'YLabel')),pltIdx);
% Add colorbar to plot 14 (2nd row 3 rd plot)
c_pos = [0.485 0.7 0.01 0.15];
colorbar(H_plots(14),'Position',c_pos,'Limits',cLims);

% Add figure annotation
oName = sprintf('sdfsHeatmap_%s_%s_%s.pdf',useFilter.outcome,useFilter.epoch,num2str(useFilter.pval*100,'Pval_%02d'));
figTitle = 'Heatmap of normalized SDFs for all units';
filterTitle = sprintf('Spike count corr. filtered for [%s, %s, and pval <= %s]',...
    useFilter.outcome,useFilter.epoch,num2str(useFilter.pval,'%0.02f'));
strs = {oName,figTitle,filterTitle};

annotation('textbox','Position',[0.01 0.96 0.98 0.03],'String',strs{1},...
    'FontSize',13,'FontWeight','bold','LineStyle','none','Interpreter','none',...
    'Color',[0.5 0.5 0.5],'FontAngle','italic');
annotation('textbox','Position',[0.31 0.97 0.98 0.03],'String',strs{2},...
    'FontSize',15,'FontWeight','bold','LineStyle','none','Interpreter','none',...
    'Color','k');
annotation('textbox','Position',[0.58 0.96 0.98 0.03],'String',strs{3},...
    'FontSize',13,'FontWeight','bold','LineStyle','none','Interpreter','none',...
    'Color',[0.5 0.5 0.5]);
% save pdfFile:
saveFigPdf(fullfile(outputPdfDir,oName));
delete(H_Figure);

end

%% Other functions
function [] = showImage(H_axis,sdfImg,timeMs,epoch,alignedEvent,unitArea,cLims)
    axes(H_axis);
    % cLims = [-1 1];
    % show image
    imagesc(sdfImg,cLims);
    % row 1 = top of image, ie the SDF of sdfImg(1,:)
    set(gca,'YDir','reverse');
    xTickLabel = unique([0:-200:min(timeMs) 0:200:max(timeMs)]);
    if ~xTickLabel(1) % first label is zero
        xTickLabel = xTickLabel(1:end-1);
    else
        xTickLabel = xTickLabel(2:end-1);
    end
    % location of ticks on the image to correspond to tick labels
    xTickLoc = arrayfun(@(x) find(timeMs==x),xTickLabel);
    % Tick location for 0 ms (align Event time)
    x0Loc = xTickLoc(xTickLabel==0);
    % add zero line
    yMax = size(sdfImg,1);
    line([x0Loc x0Loc],[0 yMax+1],'Color','k');
    % Labels etc
    title(epoch);
    h = get(gca,'YLabel');
    set(h,'String',{sprintf('%s ',unitArea),sprintf('[%d] ',size(sdfImg,1))});
    set(h,'Rotation',0,'FontWeight','bold','HorizontalAlignment','right');
    set(gca,'YTickLabel',{});
    set(get(gca,'XLabel'),'String',sprintf('Time from %s (ms)',alignedEvent)) ;
    set(gca,'XTick', xTickLoc);
    set(gca,'XTickLabel', xTickLabel);
    set(gca,'XMinorGrid','on', 'GridColor','k');
end

function [sigNonSigUnits] = getUnitNums(unitsTbl,epoch,outcome,pval)
    idx = ismember(unitsTbl.filter_Epoch,epoch) & ismember(unitsTbl.filter_Outcome,outcome) & unitsTbl.filter_Pval == pval;
    sigNonSigUnits = unitsTbl(idx,{'filter_IsRscSignificant','SEF','FEF','SC'}) ;
end

function [satSdfsTbl] = getNormalizedSdf(satSdfDir,unitNum)
    fn = fullfile(satSdfDir,num2str(unitNum,'Unit_%03d.mat'));
    sdfs = load(fn,'sdfs');
    sdfs = sdfs.sdfs;
    % sdfs.satCondition
    sdfs.satCondition = regexprep(sdfs.condition,{'Correct','Error.*'},{'',''});
    % sdfs.outcome
    sdfs.outcome = regexprep(sdfs.condition,{'Fast','Accurate'},{'',''});
    % timestamps for SDFs
    tsV = sdfs.Visual_timeMs{1};
    %% New method using Baseline Win = [-400 -200]
    % do normalization for trail sdfs for all outcomes by SatCondition
    % Compute mean SDFs for all conditions
    satConds = unique(sdfs.satCondition);
    epochs = {'Visual','PostSaccade','PostReward'};
    satSdfsTbl = table();
    for sc = 1:numel(satConds)
        tempSatSdfTbl = table();
        satCondition = satConds{sc};
        tempSatSdfTbl.unitNum = unitNum;
        tempSatSdfTbl.satCondition = {satCondition};
        idxSatCond = ismember(sdfs.satCondition,satCondition);
        tempSatSdfTbl.nTrials = sum(sdfs.Visual_nTrials(idxSatCond));
        for ep = 1:numel(epochs)
            epoch = epochs{ep};
            meanSatSdf = mean(cat(1,sdfs.([epoch '_sdfByTrial']){idxSatCond}),1);
            tempSatSdfTbl.([epoch 'AlignedEvent']) = sdfs.([epoch '_alignedEvent'])(1);
            tempSatSdfTbl.([epoch 'TimeMs']) = {sdfs.([epoch '_timeMs']){1}}; %#ok<CCAT1>
            tempSatSdfTbl.([epoch 'SatSdfMean']) = {meanSatSdf};
        end
        satSdfsTbl = [satSdfsTbl;tempSatSdfTbl]; %#ok<*AGROW>
    end
    %% Normalize SDFs by computing as follows
    % (1) frBl min(mean fr in baselineWin for Visual epoch for SAT)
    % (2) fxMax max(max(mean fr for conds SAT [F/A] by [V,PS,PR]))
    nRows = size(satSdfsTbl,1);
    baselineWin = [-400 -200];
    blIdx = find(tsV==baselineWin(1)):find(tsV==baselineWin(2));
    % get mean SDFs by SAT for 3 epochs
    frV = satSdfsTbl.VisualSatSdfMean;
    frPs = satSdfsTbl.PostSaccadeSatSdfMean;
    frPr = satSdfsTbl.PostRewardSatSdfMean;
    % (1) frBl min(mean fr in baselineWin for Visual epoch for SAT)
    frBl = min(cellfun(@(x) mean(x(blIdx)),frV));
    frBlStd = mean(cellfun(@(x) std(x(blIdx)),frV));
    % (2) fxMax max(max(mean fr for conds SAT [F/A] by [V,PS,PR]))
    frMax = max(cellfun(@max,[frV;frPs;frPr]));
    % Add baseline Fr and max Fr to output table
    satSdfsTbl.frBaseline = repmat(frBl,nRows,1);
    satSdfsTbl.frBaselineStd = repmat(frBlStd,nRows,1);
    satSdfsTbl.frMaxAllConditionsEpochs = repmat(frMax,nRows,1);
    % Compute and add Normalized sdf for SAT: normalization =
    % (frVector-BaselineFr)/(maxFr-BaselineFr)
    % Visual
    satSdfsTbl.VisualSatSdfNormalized = cellfun(@(x)(x-frBl)./(frMax-frBl),frV,'UniformOutput',false);
    % PostSaccade
    satSdfsTbl.PostSaccadeSatSdfNormalized = cellfun(@(x)(x-frBl)./(frMax-frBl),frPs,'UniformOutput',false);
    %PostReward
    satSdfsTbl.PostRewardSatSdfNormalized = cellfun(@(x)(x-frBl)./(frMax-frBl),frPr,'UniformOutput',false);

end

% Get figure template
function [H_plots,H_Figure] = getPlotHandles()
    % For Heatmap of aggregated SDFs for each area
    % 6 rows {3-Areas for Fast, 3-areas for Accurate}
    % 6 cols {3-epochs for significant, 3-epochs for non-significant}
    H_Figure = newFigure();
    nRows = 6;
    startOffset = 0.06;
    pltH = 0.125;
    pltWFor3Cols = 0.4;
    gutterx = 0.005;
    guttery = 0.010;
    % conditions SDFs have different lengths of times (x axis).
    % make widths proportional such that time unit ticks are of equal length
    % visual:[-600 400] , postSaccade:[-200 600], postReward:[-100 400]
    pltWProp = [1001 801 501];
    % partition total plot width to the proportion above
    pltWs = pltWProp.*(pltWFor3Cols/sum(pltWProp));
    offsetsX = [0 pltWs(1)+gutterx sum(pltWs(1:2))+gutterx*2] + startOffset;
    offsetsY = 0.92-pltH:-(pltH+guttery):guttery; % for 6 row-starts
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
                pos(2) = offsetsY(ro) - guttery*6;
            end
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount));
        end
    end
end
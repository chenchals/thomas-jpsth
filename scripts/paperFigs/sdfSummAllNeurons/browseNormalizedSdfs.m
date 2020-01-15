% Uses workspace variable: satSdfsTbl computed by plotSdfsByArea.m

% get all time axis vectors
tsV = satSdfsTbl.VisualTimeMs{1}';
tsPs = satSdfsTbl.PostSaccadeTimeMs{1}';
tsPr = satSdfsTbl.PostRewardTimeMs{1}';

%%
% Make subplots for [V,PS,PR]
[H_plots,H_figure] = getPlotHandles([numel(tsV),numel(tsPs),numel(tsPr)]);
% use normalized or just mean sdfs?
useSuffix = 'SatSdfNormalized';
v = 'Visual';
ps = 'PostSaccade';
pr = 'PostReward';
fClr = [0 0.7 0];
aClr = [1 0 0];
unitNums = unique(satSdfsTbl.unitNum,'stable');
for un = 1:numel(unitNums)
    unitNum = unitNums(un);
    sdfs = satSdfsTbl(satSdfsTbl.unitNum==unitNum,:);
    fastIdx = find(ismember(sdfs.satCondition,'Fast'));
    accuIdx = find(ismember(sdfs.satCondition,'Accurate'));
    titleStr = sprintf('%s Unit# %d',sdfs.unitArea{1},unitNum);
    if sdfs.isRscSignificant
        titleStr = [titleStr '   [Rsc Significant]']; %#ok<*AGROW>
    else
        titleStr = [titleStr '   [Rsc *not* Significant]']; %#ok<*AGROW>
    end
    clrs = [fClr;aClr];
    % get SDFs for all epochs
    frV = cell2mat(sdfs.([v useSuffix])([fastIdx,accuIdx]))';
    frPs = cell2mat(sdfs.([ps useSuffix])([fastIdx,accuIdx]))';
    frPr = cell2mat(sdfs.([pr useSuffix])([fastIdx,accuIdx]))';
    % update plots
    updatePlot(H_plots(1),tsV,frV,v,'CueOn',clrs);
    updatePlot(H_plots(2),tsPs,frPs,ps,'SaccadePrimary',clrs);
    updatePlot(H_plots(3),tsPr,frPr,pr,'Reward',clrs);
    % update yAxis
    yMinMax = round(minmax([frV(:);frPs(:);frPr(:)]'),2);
    yMinMax = yMinMax + [-0.05 0.05];
    set(H_plots,'YLim',yMinMax);
    set(get(H_plots(1),'YLabel'),'String',useSuffix);
    set(H_plots(2:3),'YTickLabel',{});
    % add annotation
    h_a = annotation('textbox','String',titleStr,'FontSize',18,'FontWeight','bold',...
        'Position',[0.3 0.98 0.5 0.04],'LineStyle','none');
    pause
    % do nto delete annotation if it is the last unit
    if un < numel(unitNums)
      delete(h_a)
    end
end

%%
function [H_axis] = updatePlot(H_axis,ts,fr,epoch,alignedOn,clrs)
   axes(H_axis);
   h=plot(ts,fr);
   h(1).Color = clrs(1,:);
   h(2).Color = clrs(2,:);
   grid('on')
   xlabel(sprintf('Time from %s [ms]',alignedOn));
   title(epoch);
end

function [H_plots,H_figure] = getPlotHandles(pltWProp)
    H_figure = newFigure();
    pos = [0.05 0.65 0.92 0.3];
    set(H_figure,'Position',pos);
    startOffset = 0.05;
    pltH = 0.75;
    %total plot width to the proportion above
    pltWFor3Cols = 0.9;
    gutterx = 0.005;
    % conditions SDFs have different lengths of times (x axis).
    % make widths proportional such that time unit ticks are of equal length
    % visual:[-600 400] , postSaccade:[-200 600], postReward:[-100 400]
    % partition total plot width to the proportion above
    pltWs = pltWProp.*(pltWFor3Cols/sum(pltWProp));
    offsetsX = [0 pltWs(1)+ gutterx sum(pltWs(1:2))+gutterx*2] + startOffset;
    offsetsY = 0.9-pltH; % for 1 row
    pltCount = 0;
    for col = 1:numel(offsetsX)
        pos(1) = offsetsX(col);
        pos(2) = offsetsY;
        pos(3:4) = [pltWs(col) pltH];
        pltCount = pltCount + 1;
        H_plots(pltCount) = axes('parent',H_figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount));
    end

end

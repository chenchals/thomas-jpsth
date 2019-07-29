function [plotHandles] = jpsthFigTemplate(datafile)
%function [plotHandles] = jpsthFigTemplate(datafile,targetLocs,conditions,alignedOn)
% 4 jpsths figure
%|             | fastCorrect | accurateCorrect|
%|-------------|-------------|----------------|
%| targetLoc#1 |     1,1     |      1,2       |
%|-------------|-------------|----------------|
%| targetLoc#2 |     2,1     |      2,2       |
%|---------------------------|----------------|

%datafile ='/Volumes/schalllab/Users/Chenchal/JPSTH/FEF_SC_Visual_1ms/PAIR_0077_D20130828001-RH_SEARCH_DSP12a_FEF_DSP17a_SC.mat';
% collect data and labels:
targetLocs= {'TargetInXandY' 'TargetInXnotY'};
conditions = {'AccurateCorrect' 'FastCorrect'};
alignedOn = 'CueOn';

temp = load(datafile);
nJpsths = 0;
for c=1:2
    for r=1:2
        nJpsths = nJpsths + 1;
        jpsthLabels{nJpsths} = join({targetLocs{r},conditions{c}},'-');
        jpsths{nJpsths} = temp.(targetLocs{r}).(conditions{c})(alignedOn,:);
    end
end
jpsths = jpsths(end:-1:1);
jpsthLabels = jpsthLabels(end:-1:1);
cellPairInfo = temp.cellPairInfo;
singletonLocs = temp.singletonLocs;
clear temp;

%% Compute positions of JPSTH image on normalized axes, in arbitarary units
coinsXoffset = 3;
boxcarFilt = 5; 
allPlotsShiftX = 20;
allPlotsShiftY = -10;

%% JPSTH position - compute other plot pos based on this position
jpsthWH = 100;
% from bottom col1, row1
jpsthPos{1} = [40+allPlotsShiftX 60+allPlotsShiftY jpsthWH jpsthWH];
% from bottom col1, row2
jpsthPos{2} = [jpsthPos{1}(1) jpsthPos{1}(2)+2*jpsthWH jpsthWH jpsthWH];
% from bottom col2, row1
jpsthPos{3} = [jpsthPos{1}(1)+3*jpsthWH jpsthPos{1}(2) jpsthWH jpsthWH];
% from bottom col2, row2
jpsthPos{4} = [jpsthPos{1}(1)+3*jpsthWH jpsthPos{1}(2)+2*jpsthWH jpsthWH jpsthWH];

%% Coincidence Histogram plots 
% rotated and aligned to bottom of jpsth and to the right. The rorated plot
% is zoomed to sqrt(2) so that the physical length of the diagonal
% corresponds to the x-axis lims.
coinsPos = cellfun(@(x) [x(1)+coinsXoffset+x(3) x(2) x(3) x(4)],jpsthPos,'UniformOutput',false);

%% XCell PSTH plots
psthHeight = 30;
psthOffset = 3;
% below jpsth
xPsthPos1 = cellfun(@(x) [x(1) x(2)-psthHeight-psthOffset x(3) psthHeight],jpsthPos,'UniformOutput',false);
% below coincidence hist
xPsthPos2 = cellfun(@(x) [x(1)+x(3)+coinsXoffset x(2)-psthHeight-psthOffset x(3) psthHeight],jpsthPos,'UniformOutput',false);

%% YCell PSTH plots
% Left of jpsth
yPsthPos = cellfun(@(x) [x(1)-psthHeight-psthOffset x(2) psthHeight x(4)],jpsthPos,'UniformOutput',false);

%% CrossCorrelation Histogram plots
% rotated and centered at right-top vertex fo coincidence histogram
% plot clipped to half of x-axis data ie if JPSTH is [-300 300] --> then
% crosscorr is clipped to [-150 150], as the polt is not zoomed to sqrt(2)
% after rotation, so the physical length of the diagonal corresponds to half the x-axis 
xcorrNudge = hypot(coinsXoffset,psthOffset)/2 ;
xcorrOffset = 5+xcorrNudge;
xcorrsPos = cellfun(@(x) [xcorrOffset+xcorrNudge+x(1)+x(4)*1.5 xcorrOffset-xcorrNudge+x(2)+x(3)*0.5 x(3) x(4)],jpsthPos,'UniformOutput',false);

%% Draw new Figure
delete(allchild(0))
figureName='JPSTH Analyses';

fontName = 'Arial';
fontSize = 8;
lineWidth = 0.1;
screenWH=get(0,'ScreenSize');%pixels
screenPixPerInch=get(0,'ScreenPixelsPerInch');%pixels
screenSize = [screenWH(3) screenWH(4)]./screenPixPerInch;
% figure in pixels
fW = 1600; 
fH = 1200;
figurePos = [20 20 fW fH];

set(0,'units','pixels')
set(0,'defaulttextfontsize',fontSize,...
    'defaultaxesfontsize',fontSize,...
    'defaulttextfontname',fontName,...
    'defaultaxeslinewidth',lineWidth)
%Main figure window
H_Figure=figure('Position',figurePos,...
    'PaperOrientation', 'Landscape',...
    'NumberTitle','off',...
    'Units','pixels',...
    'Menu','none',...
    'Name',figureName,...
    'Color',[1 1 1],...
    'Tag','Figure');

set(H_Figure,'Units','normalized');
figPos = get(H_Figure,'Position');
% Scale the arbitarary units of position on to normalized axes
% assume 0-1 x axis correponds to 1000 Units
sW = 1000/fW;
sH = 1000/fH;
fx_getPos = @(p) [p(1)*sW p(2)*sH p(3)*sW p(4)*sH].*1/400;
grayCol = [0.6 0.6 0.6];
fastCol = [0 1 0];
accurateCol = [1 0 0];

% JPSTHs
jpsthMinMax = cell2mat(cellfun(@(x) minmax(x.normalizedJpsth{1}(:)'),jpsths,'UniformOutput',false));
jpsthMinMax = minmax(jpsthMinMax(:)');
jpsthMinMax = round(jpsthMinMax.*10)./10;

% Coincidences
coinsMinMax = cell2mat(cellfun(@(x) minmax(x.coincidenceHist{1}(:,2))',jpsths,'UniformOutput',false));
coinsMinMax = minmax(coinsMinMax(:)');
coinsAbsMax = max(abs(coinsMinMax));
coinsYlim = [-coinsAbsMax coinsAbsMax].*sqrt(2);
coinsXlim = minmax(jpsths{1}.coincidenceHist{1}(:,1)');

% XCorrelations
crossScaleY = 2;
% Rotating by 45 deg, reduces the physical viewable x axis to half when aligned on non-rotated plotbox
% so change xLim to half the data centerd at 0
crossScaleX = 0.5; 
crossMinMax = cell2mat(cellfun(@(x) minmax(x.xCorrHist{1}(:,2)'),jpsths,'UniformOutput',false));
crossMinMax = minmax(crossMinMax(:)');
crossAbsMax = max(abs(crossMinMax));
crossYlim = [-crossAbsMax crossAbsMax].*crossScaleY;
crossXlim = minmax(jpsths{1}.xCorrHist{1}(:,1)').*crossScaleX;

frFactor = 1000/jpsths{1}.binWidth;
% XCell PSTH
xCellFrMax = cellfun(@(x) max(x.xPsth{1}),jpsths,'UniformOutput',false)';
% YCell PSTH
yCellFrMax = cellfun(@(x) max(x.yPsth{1}),jpsths,'UniformOutput',false)';
% All PSTH xlims
psthLim = minmax(jpsths{1}.xPsthBins{1});

for ii = 1:4
    label = jpsthLabels{ii};
    if contains(label,'Fast')
        histColor = fastCol;
        coinsColor = fastCol;
    elseif contains(label,'Accurate')
        histColor = accurateCol;
        coinsColor = accurateCol;
    else
        histColor = grayCol;
        coinsColor = grayCol;
    end
    
    %% JPSTH
    pos = fx_getPos(jpsthPos{ii});
    H_jpsth(ii) = axes(H_Figure,'Position',pos,'Box','on','Tag','H_Figure');
    imagesc(flipud(jpsths{ii}.normalizedJpsth{1}),jpsthMinMax);
    set(H_jpsth(ii),'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
    %% Coincidence Histogram
    pos = fx_getPos(coinsPos{ii});
    H_coins(ii) = axes(H_Figure,'Position',pos,'Box','off','Tag',['H_coins' num2str(ii,'_%d')]);
    set(bar(jpsths{ii}.coincidenceHist{1}(:,1),smooth(jpsths{ii}.coincidenceHist{1}(:,2)),boxcarFilt),...
        'Facecolor',coinsColor,'FaceAlpha', 0.2,'edgecolor','none');
    line(coinsXlim,[0 0],'color','k')
    set(H_coins(ii),'XLim',coinsXlim,'YLim',coinsYlim,'Box', 'off');
    axis off
    camzoom(sqrt(2));
    camorbit(-45,0);
    %% XCell PSTH
    % xpsth1
    pos = fx_getPos(xPsthPos1{ii});
    H_xpsth1(ii) = axes(H_Figure,'Position',pos,'Box','on','Tag',['H_xpsth1' num2str(ii,'_%d')]);
    set(bar(jpsths{ii}.xPsthBins{1},smooth(jpsths{ii}.xPsth{1}),boxcarFilt),...
        'Facecolor',histColor,'FaceAlpha', 0.2,'edgecolor','none');
    set(H_xpsth1(ii),'YDir','Reverse'); 
    set(H_xpsth1(ii),'XLim',psthLim,'YLim',[0 xCellFrMax{ii}]);
    % xpsth2
    pos = fx_getPos(xPsthPos2{ii});
    H_xpsth2(ii) = axes(H_Figure,'Position',pos,'Box','on','Tag',['H_xpsth2' num2str(ii,'_%d')]);
    set(bar(jpsths{ii}.xPsthBins{1},smooth(jpsths{ii}.xPsth{1}),boxcarFilt),...
        'Facecolor',histColor,'FaceAlpha', 0.2,'edgecolor','none');
    set(H_xpsth2(ii),'YDir','Reverse');
    set(H_xpsth2(ii),'XLim',psthLim,'YLim',[0 xCellFrMax{ii}],'YAxisLocation','right');
    %% YCell PSTH
    pos = fx_getPos(yPsthPos{ii});
    H_ypsth1(ii) = axes(H_Figure,'Position',pos,'Box','on','Tag',['H_jpsth1' num2str(ii,'_%d')]);
    set(bar(jpsths{ii}.yPsthBins{1},smooth(jpsths{ii}.yPsth{1}),boxcarFilt),...
        'Facecolor',histColor,'FaceAlpha', 0.2,'edgecolor','none');
    set(H_ypsth1(ii),'YDir','Reverse'); 
    set(H_ypsth1(ii),'XLim',psthLim,'YLim',[0 yCellFrMax{ii}]);
    view([90 -90])
    %% Cross Correlation Histogram
    pos = fx_getPos(jpsthPos{ii});
    H_xcorr(ii) = axes(H_Figure,'Position',pos,'Box','off','Tag',['H_xcorr' num2str(ii,'_%d')]);
    set(bar(jpsths{ii}.xCorrHist{1}(:,1),smooth(jpsths{ii}.xCorrHist{1}(:,2)),boxcarFilt),...
        'Facecolor',grayCol,'edgecolor','none');
    line(crossXlim,[0 0],'color','k')
    line([0 0],[0 crossAbsMax],'color','k')
    set(H_xcorr(ii),'XLim',crossXlim,'YLim',crossYlim,'Box', 'off');
        axis off
    %camzoom(sqrt(2));
    camorbit(45,0);
    set(H_xcorr(ii),'Position',fx_getPos(xcorrsPos{ii}));     
   drawnow
end

% Annotate:
% add plot labels
addJpsthLabel(jpsthLabels,H_ypsth1);
annotateYCellPsth(H_ypsth1,frFactor,cellPairInfo);
annotateXCellPsth(H_xpsth1,frFactor,cellPairInfo);
annotateXCellPsth(H_xpsth2,frFactor,cellPairInfo);


vars = who;
vars= vars(contains(vars,'H_'));
for ii = 1:numel(vars)  
    if ishandle(eval(vars{ii}))
        plotHandles.(vars{ii}) = eval(vars{ii});
    end
end
end

function addJpsthLabel(labels,axesHandles)
  for ii = 1:numel(axesHandles)
    h=axesHandles(ii);
    text(h,max(get(h,'XLim'))*1.4,max(get(h,'YLim'))*1.1,labels{ii},...
        'FontWeight','bold','FontSize',14,'fontAngle','Italic');
  end
end

function annotateYCellPsth(axesHandles,frFactor,cellPairInfo)
  for ii = 1:numel(axesHandles)
    xLims = get(axesHandles(ii),'XLim');
    yLims = get(axesHandles(ii),'YLim').*1.01;
    xTicks = min(xLims):100:max(xLims);
    yTickLabel = num2str(round(max(yLims)*frFactor),'%d spk/s');
    set(axesHandles(ii),'XTick',xTicks(2:end-1),'XTickLabelRotation',90,'XGrid','on');    
    set(axesHandles(ii),'YTick',[]);
    text(axesHandles(ii),min(xLims),max(yLims)*1.1,yTickLabel,'Rotation',135,...
        'HorizontalAlignment','center','VerticalAlignment','baseline')
  end
end

function annotateXCellPsth(axesHandles,frFactor,cellPairInfo)
for ii = 1:numel(axesHandles)
    xLims = get(axesHandles(ii),'XLim');
    yLims = get(axesHandles(ii),'YLim').*1.01;
    xTicks = min(xLims):100:max(xLims);
    yTickLabel = num2str(round(max(yLims)*frFactor),'%d spk/s');
    set(axesHandles(ii),'XTick',xTicks(2:end-1),'XTickLabelRotation',90,'XGrid','on');
    set(axesHandles(ii),'YTick',[]);
    if strcmpi(get(axesHandles(ii),'YAxisLocation'),'right')
         text(axesHandles(ii),max(xLims),max(yLims)*1.1,yTickLabel,'Rotation',45,...
        'HorizontalAlignment','center','VerticalAlignment','baseline')
    else   
         text(axesHandles(ii),min(xLims),max(yLims)*1.1,yTickLabel,'Rotation',135,...
        'HorizontalAlignment','center','VerticalAlignment','baseline')
    end
end
end

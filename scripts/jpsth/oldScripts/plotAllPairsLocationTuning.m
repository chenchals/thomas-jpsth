
% Plot tuning curves for all pairs of cells 
trialTypesDb = load('/Volumes/schalllab/Users/Chenchal/JPSTH/TrialTypesDB.mat');
singletonLocsBySession = trialTypesDb.TrialTypesDB(:,{'session','SingletonLoc'});
figDir = '/Volumes/schalllab/Users/Chenchal/JPSTH/FigsBinWidth_5';%
tuneDir = fullfile(figDir,'tuning');
if ~exist(tuneDir,'dir')
    mkdir(tuneDir);
end
plotlineStyle.FEF = 'k--';
plotlineStyle.SEF = 'k:';
plotlineStyle.SC = 'k-';

files = dir(fullfile(figDir,'PAIR*.mat'));
fns = {files.name}';
files = strcat({files.folder}',filesep,{files.name}');
conditions = {'AccurateCorrect'; 'FastCorrect'};
cueMeanWindow = [50 200];
saccMeanWindow = [-100 100];
nPairs = numel(files);

H_Plots = createFigureTemplate();
nextPlot = 0;
fileNum = 0;
for cp = 1:nPairs
    pairData = load(files{cp},conditions{:},'cellPairInfo');
    cellPairInfo = pairData.cellPairInfo;
    sessionName = char(regexprep(cellPairInfo.datafile,'-.*$',''));
    locsByTrialNo = singletonLocsBySession.SingletonLoc{strcmpi(singletonLocsBySession.session,sessionName)};
   
    infoText{1} = {'CellId';pairData.cellPairInfo.X_cellIdInFile{1};pairData.cellPairInfo.Y_cellIdInFile{1}};
    infoText{2} = {'depth';num2str(pairData.cellPairInfo.X_depth);num2str(pairData.cellPairInfo.Y_depth)};
    infoText{3} = {'area';pairData.cellPairInfo.X_area{1};pairData.cellPairInfo.Y_area{1}};   
    infoText{4} = {'RF';['[' num2str(pairData.cellPairInfo.X_RF{1}) ']'];['[' num2str(pairData.cellPairInfo.Y_RF{1}) ']']};
    infoText{5} = {'MF';['[' num2str(pairData.cellPairInfo.X_MF{1}) ']'];['[' num2str(pairData.cellPairInfo.Y_MF{1}) ']']};
    
    lineStyles = {plotlineStyle.(pairData.cellPairInfo.X_area{1}); plotlineStyle.(pairData.cellPairInfo.Y_area{1})};
    
    % Aggregate all corrects
    binWidth = pairData.AccurateCorrect{'CueOn','binWidth'};
    correctTrialNos = cell2mat([pairData.AccurateCorrect{'CueOn','trialNosByCondition'};pairData.FastCorrect{'CueOn','trialNosByCondition'}]);
    correctSingletonLocs = locsByTrialNo(correctTrialNos);

    %% Aligned on CueOn
    cueOnTimeBins = cell2mat(pairData.AccurateCorrect{'CueOn','xPsthBins'})';
    xCueOnSpkCount = cell2mat([pairData.AccurateCorrect{'CueOn','xSpikeCounts'};pairData.FastCorrect{'CueOn','xSpikeCounts'}]);
    yCueOnSpkCount = cell2mat([pairData.AccurateCorrect{'CueOn','ySpikeCounts'};pairData.FastCorrect{'CueOn','ySpikeCounts'}]);
    % CueOn SpikeCounts by location over averaging window
    cueIdx = cueOnTimeBins>=cueMeanWindow(1) & cueOnTimeBins<=cueMeanWindow(2);
    xCueOnByLocs = arrayfun(@(x) xCueOnSpkCount(correctSingletonLocs==x-1,cueIdx),1:8,'UniformOutput',false);
    yCueOnByLocs = arrayfun(@(x) yCueOnSpkCount(correctSingletonLocs==x-1,cueIdx),1:8,'UniformOutput',false);
    
    xCueOnMeanFrByLoc = cellfun(@(x) mean(x(:))*1000/binWidth,xCueOnByLocs);
    yCueOnMeanFrByLoc = cellfun(@(x) mean(x(:))*1000/binWidth,yCueOnByLocs);
    %% Aligned on SaccadePrimary
    saccadePrimaryTimeBins = cell2mat(pairData.AccurateCorrect{'SaccadePrimary','xPsthBins'})';
    xSaccadePrimarySpkCount = cell2mat([pairData.AccurateCorrect{'SaccadePrimary','xSpikeCounts'};pairData.FastCorrect{'SaccadePrimary','xSpikeCounts'}]);
    ySaccadePrimarySpkCount = cell2mat([pairData.AccurateCorrect{'SaccadePrimary','ySpikeCounts'};pairData.FastCorrect{'SaccadePrimary','ySpikeCounts'}]);
    % SaccadePrimary SpikeCounts by location
    saccIdx = saccadePrimaryTimeBins>=saccMeanWindow(1) & saccadePrimaryTimeBins<=saccMeanWindow(2);
    xSaccPrimaryByLocs = arrayfun(@(x) xSaccadePrimarySpkCount(correctSingletonLocs==x-1,saccIdx),1:8,'UniformOutput',false);
    ySaccPrimaryByLocs = arrayfun(@(x) ySaccadePrimarySpkCount(correctSingletonLocs==x-1,saccIdx),1:8,'UniformOutput',false);
 
    xSaccPrimaryMeanFrByLoc = cellfun(@(x) mean(x(:))*1000/binWidth,xSaccPrimaryByLocs);
    ySaccPrimaryMeanFrByLoc = cellfun(@(x) mean(x(:))*1000/binWidth,ySaccPrimaryByLocs);
    
    %% Plot tuning
    nextPlot = nextPlot + 1;
    axes(H_Plots(nextPlot));
    plot([0:7],xSaccPrimaryMeanFrByLoc,lineStyles{1},'LineWidth',2);
    hold on
    plot([0:7],ySaccPrimaryMeanFrByLoc,lineStyles{2},'LineWidth',2);
    hold off
    %plot([0:7],[xSaccPrimaryMeanFrByLoc;ySaccPrimaryMeanFrByLoc]);
    xlabel('Singleton Location Index', 'FontSize', 12, 'FontWeight','bold');
    ylabel(['Mean FR in [' num2str(saccMeanWindow) '] ms window (Saccade Primary)'], 'FontSize', 12, 'FontWeight','bold');
    for t = 1:5
       text(t,max(ylim),infoText{:,t},'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');
    end
    text(range(xlim)/2,max(ylim),strcat(cellPairInfo.Pair_UID, ' : ', sessionName),'FontSize', 12, 'FontWeight','bold','Interpreter','none','VerticalAlignment','bottom','HorizontalAlignment','center');
   drawnow;

    if mod(cp,12) == 0
        fileNum = fileNum +1;
        fn = fullfile(tuneDir,num2str(fileNum,'Tuning_SaccadePrimary_aligned_%02d'));
        saveFigAs(fn);
        nextPlot = 0;
        arrayfun(@cla,H_Plots);
    end
     
end

function saveFigAs(fn)
    % currUnits = get(gcf,'Units');
    % currPosition = get(gcf,'Position');
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    fprintf('Saving figure to: %s\n',fn); 
    print(fn,'-dpdf','-painters')
    drawnow
end

function [H_Plots] = createFigureTemplate()
    %create figure template
    delete(allchild(0))
    figureName='Burst Analyses';

    fontName = 'Arial';
    fontSize = 11;
    lineWidth = 0.05;
    screenWH=get(0,'ScreenSize');%pixels
    margin = 40;
    figureWidth = screenWH(3)-2*margin;
    figureHeight = screenWH(4)-2*margin;
    figurePos = [margin margin figureWidth figureHeight];

    set(0,'units','pixels')
    set(0,'defaulttextfontsize',fontSize,...
        'defaultaxesfontsize',fontSize,...
        'defaulttextfontname',fontName,...
        'defaultaxeslinewidth',lineWidth)

    %Main figure window
    H_Figure=figure('Position',figurePos,...
        'PaperOrientation', 'Landscape',...
        'NumberTitle','off',...
        'Menu','none',...
        'Units','normalized',...
        'Name',figureName,...
        'Color',[1 1 1],...
        'Tag','H_Figure');
    % OutlinePos=[0.005 0.005 0.99 0.99];
    % H_Outline=axes('parent',H_Figure,...
    %                'position',OutlinePos,...
    %                'box','on',...
    %                'xtick',[],...
    %                'ytick',[],...
    %                'xcolor',[0 0 0],...
    %                'ycolor',[0 0 0],...
    %                'Tag','H_Outline');

    xGutter = 0.02;
    yGutter = 0.04;
    nCols = 4;
    plotW = (1/nCols) - 1.5 * xGutter;
    nRows = 3;
    plotH = (1/nRows) - 1.5 * yGutter;
    plotX = ((0:nCols-1).*(1/nCols)) + xGutter;
    plotX = repmat(plotX,nRows,1);
    plotY = ((0:nRows-1).*(1/nRows)) + yGutter;
    plotY = repmat(plotY(:),1,nCols);
    
    H1 = arrayfun(@(x) axes('Parent',H_Figure,...
        'Position',[plotX(3,x) plotY(3,x) plotW plotH],...
        'Box','on',...
        'Tag',['H_' num2str(x)]),...
        (1:4),'UniformOutput',false);
    H2 = arrayfun(@(x) axes('Parent',H_Figure,...
        'Position',[plotX(2,x) plotY(2,x) plotW plotH],...
        'Box','on',...
        'Tag',['H_' num2str(x+4)]),...
        (1:4),'UniformOutput',false);
    H3 = arrayfun(@(x) axes('Parent',H_Figure,...
        'Position',[plotX(1,x) plotY(1,x) plotW plotH],...
        'Box','on',...
        'Tag',['H_' num2str(x+8)]),...
        (1:4),'UniformOutput',false);
    H_Plots = [H1{:} H2{:} H3{:}];
    
end
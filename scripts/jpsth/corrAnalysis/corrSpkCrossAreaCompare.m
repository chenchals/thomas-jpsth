% Compare cross area and within area spike count correlations
% Use static windows
fn = 'dataProcessed/analysis/spkCorr/spkCorrAllPairsStatic.mat';
spkCorr = load(fn);
outPdfFn = 'dataProcessed/analysis/spkCorr/Summary_SpkCorr';

%% 
filterData = 1;
filterIdx = 2;

% group stats for the following areas and align windows
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
% Prune to columns of interest for each pairArea
useCols = {
    'X_monkey'
    'pairAreas'
    'nTrials'
    'XY_Dist'
    'condition'
    'alignedName'
    'alignedEvent'
    'rho_pval_win'
    'rhoRaw'
    'pvalRaw'
    'rhoZTrial'
    'pvalZTrial'
    'signifRaw_05'
    'signifRaw_01'
    'signifZTrial_05'
    'signifZTrial_01'
    };
% compareConditions = {
%     'AccurateCorrect' 'FastCorrect'
%     'AccurateErrorChoice' 'FastErrorChoice'
%     'AccurateErrorTiming' 'FastErrorTiming'
%     };
% aggregate required fields into a single table
spkCorrAll = table();
for pa =1:numel(pairAreas)
    pairArea = pairAreas{pa};
    spkCorrAll = [spkCorrAll; spkCorr.(pairArea)(:,useCols)];
end
% recode X_monkey
spkCorrAll.monkey = spkCorrAll.X_monkey;
spkCorrAll.X_monkey = [];
% use absolute values
spkCorrAll.rhoRaw = abs(spkCorrAll.rhoRaw);
spkCorrAll.rhoZTrial = abs(spkCorrAll.rhoZTrial);
spkCorrAll.baseCondition = regexprep(spkCorrAll.condition,'(Accurate)|(Fast)','');
spkCorrAll.trialType = regexprep(spkCorrAll.condition,'(Correct)|(Error.*)','');
% Recode distances as pairAreas
spkCorrAll.newPairAreas = spkCorrAll.pairAreas;
spkCorrAll.newXY_Dist = spkCorrAll.XY_Dist;
idx = spkCorrAll.XY_Dist==0;
spkCorrAll.newPairAreas(idx) = strcat(spkCorrAll.pairAreas(idx),'(0mm)');
idx = spkCorrAll.XY_Dist>0 & spkCorrAll.XY_Dist<2;
spkCorrAll.newPairAreas(idx) = strcat(spkCorrAll.pairAreas(idx),'(<2mm)');

% use newPairAreas
newPairAreas =  {'SEF-SEF(0mm)','SEF-SEF(<2mm)','SEF-SEF','SEF-FEF','SEF-SC'};

% Filter data?
filterTxt = {};
outFileSuffix = '';
if filterData
    filters.monk = {'D','E'};    
    filt = filters.monk{filterIdx};
    filterTxt = [filterTxt, ['Monk: ' filt]];
    outFileSuffix = [outFileSuffix '_' filt];
    filteredSpkCorrAll = spkCorrAll(strcmp(spkCorrAll.monkey,filt),:);
else
    outFileSuffix = '_All';
    filteredSpkCorrAll = spkCorrAll;
end

% get stats
spkCorrAllStats = grpstats(filteredSpkCorrAll(:,{'baseCondition','trialType','newPairAreas','condition','alignedName','rhoRaw','rhoZTrial'}),...
        {'condition','baseCondition','trialType','alignedName','newPairAreas'},{'mean','std','min','max'});
spkCorrAllStats.sem_rhoRaw = (spkCorrAllStats.std_rhoRaw.^2)./sqrt(spkCorrAllStats.GroupCount);
spkCorrAllStats.sem_rhoZTrial = (spkCorrAllStats.std_rhoZTrial.^2)./sqrt(spkCorrAllStats.GroupCount);

%% Plot it
if filterData
    txtFilter = ['Filters : ' char(join(filterTxt',' ; '))];
else
    txtFilter = '';
end
pdfFn = [outPdfFn outFileSuffix '.pdf'];
plotCorrStats(pdfFn,spkCorrAllStats,newPairAreas,txtFilter);
    

%% Other functions  

function [] = plotCorrStats(pdfFn,spkCorrAllStats,newPairAreas,txtFilter)

%% Plot the stats
    baseConditions = {'Correct','ErrorChoice','ErrorTiming'};
    alignedNames = {'Baseline','Visual','PostSaccade','PostReward'};

    H_plots = getPlotHandles();
    H_plt_Idx = 0;
    for co = 1:numel(alignedNames)
        alignedName = alignedNames{co};
        for ro = 1:numel(baseConditions)
            baseCond = baseConditions{ro};
            tblRows = strcmp(spkCorrAllStats.alignedName,alignedName) & ...
                strcmp(spkCorrAllStats.baseCondition,baseCond);
            pairAreasInData = unique(spkCorrAllStats.newPairAreas(tblRows),'stable');
            pltTbl = table();
            pltTbl.datAreas = pairAreasInData;
            pltTbl.datMean = reshape(spkCorrAllStats.mean_rhoRaw(tblRows),numel(pairAreasInData),[]);
            pltTbl.datStd = reshape(spkCorrAllStats.std_rhoRaw(tblRows),numel(pairAreasInData),[]);
            pltTbl.datSem = reshape(spkCorrAllStats.sem_rhoRaw(tblRows),numel(pairAreasInData),[]);
            pltTbl.datCounts = reshape(spkCorrAllStats.GroupCount(tblRows),numel(pairAreasInData),[]);
            
            % It is possible that some pair-areas or groups may not exist
            % in data.. In these cases, add 'NaN' to data
            %if numel(pairAreasInData)~=numel(newPairAreas)               
               tempTbl = table();
               tempTbl.datAreas = newPairAreas';
               tempTbl.grpOrder = [1:size(tempTbl,1)]';
               pltTbl = outerjoin(tempTbl,pltTbl,'MergeKeys',true);
               pltTbl = sortrows(pltTbl,'grpOrder');
            %end
                        
            % plot it
            H_plt_Idx = H_plt_Idx + 1;
            plotGroupedBarsWithErrBars(H_plots(H_plt_Idx),pltTbl.datAreas,...
                pltTbl.datMean,pltTbl.datStd,pltTbl.datSem,pltTbl.datCounts);
        end
    end
    addAnnotations(pdfFn,baseConditions,alignedNames,txtFilter);
    saveFigPdf(pdfFn);

end


function addAnnotations(pdfFile,rowNames,colNames,txtFilter)
    % 
    fontSize = 18;%24
    [~,fn,ext] = fileparts(pdfFile);
    annotation('textbox',[0.15 0.97 0.05 0.03],'String',[fn ext],'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    annotation('textbox',[0.50 0.97 0.25 0.03],'String',txtFilter,'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    % conditions / alignNames
    annotation('textbox',[0.14 0.92 0.05 0.05],'String',colNames{1},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.38 0.92 0.05 0.05],'String',colNames{2},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.62 0.92 0.05 0.05],'String',colNames{3},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.86 0.92 0.05 0.05],'String',colNames{4},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    % row names    
    xOff = 0.02;
    yOff = [0.8 0.5 0.2];
    for ii = 1:numel(rowNames)
         axes('Parent',gcf,'Position',[xOff yOff(ii) 0.005 0.005],'Box','off','Visible','off')
         text(0,0,rowNames{ii},'FontSize',14,'FontWeight','bold','Rotation',90)
    end
end


function [] = plotGroupedBarsWithErrBars(H_axes,xData,yMean,yStd,ySem,nPairs)
    useSem = 1;
    axes(H_axes);
    accColor = [1 0 0];
    fastColor = [0 0.5 0];
    ylims = [0 0.50];
    %x = categorical(xData(:,1));
    x = 1:numel(xData(:,1));
    hBar = bar(x,yMean,'FaceAlpha',0.4,'BarWidth',0.90);
    set(hBar(1),'FaceColor',accColor); set(hBar(2),'FaceColor',fastColor);
    set(gca,'xticklabels',xData(:,1),'TickLabelInterpreter','none','XTickLabelRotation',20);
    set(gca,'YGrid','on')
    ylim(ylims);
    ylabel('Mean r_{sc} \pm SEM (abs.)');
    set(gca,'FontWeight','bold','FontSize',8);
    hold on
    for k1 = 1:size(yMean,2)
        ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    end
    if useSem
        yErr = ySem;
    else
        yErr = yStd;
    end
    %
    he = errorbar(ctr,yMean',yErr','YPositiveDelta',[],'Marker','o','MarkerSize',5,'LineStyle',...
        'none','LineWidth',0.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
    % add text values to plot
    [txtPairArea,txtAcc,txtFast] = arrayfun(@(x) deal(sprintf('%14s',xData{x,1}),...
        sprintf('%0.2f\\pm%0.2f %4d',yMean(x,1),yErr(x,1),nPairs(x,1)),...
        sprintf('%0.2f\\pm%0.2f %4d',yMean(x,2),yErr(x,2),nPairs(x,2))),...
        (1:size(yMean,1))','UniformOutput',false);
    txtPairArea = [sprintf('%14s','AREA');txtPairArea];
    txtAcc = [sprintf('%14s','\mu\pmSEM   n');txtAcc];
    txtFast = [sprintf('%14s','\mu\pmSEM   n');txtFast];
    
    
    y = ylims(2)-0.05;
    h_pa = text(ctr(2,2)+0.5,y,txtPairArea,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8);
    ext = get(h_pa,'Extent');
    pos = get(h_pa,'Position');
    h_txt1 = text(pos(1)+ext(3),y,txtAcc,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'color',accColor);
    ext = get(h_txt1,'Extent');
    pos = get(h_txt1,'Position');
    h_txt2 = text(pos(1)+ext(3)+0.1,y,txtFast,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'color',fastColor);
    hold off
end


function [H_plots] = getPlotHandles()
    %% for spike corr comparision among areas
    % 3 rows : correct, errorChoice, errorTiming
    % 4 cols : baseline, visual, postSaccade, postReward
    ss = [20 20 1500 990];
    margin = 10;
    FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
    %Main figure window
    H_Figure=figure('Position',FigPos,...
        'color',[1 1 1],'numbertitle','off','renderer','painters',...
        'renderermode','manual','menubar','none',...
        'Tag','H_Figure');
    orient landscape
    set(H_Figure,'units','normalized')
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);

    nRows = 3;
    nCols = 4;

    pltH = 0.25;
    pltW = 0.20;
    gutter = 0.04;
    offsetsX = 0.06:(pltW+gutter):1-gutter; % for 4 column-starts
    offsetsY = 0.95-pltH:-(pltH+gutter+0.01):gutter; % for 3 row-starts

    pltCount = 0;
    for cols = 1:nCols
        pos(1) = offsetsX(cols);
        pos(3:4) = [pltW pltH];
        for ros = 1:nRows
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
        end
    end


end
    
    

% Compare cross area and within area spike count correlations
% Use static windows
fn = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat';
%spkCorr = load(fn);
outPdfFn = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/crossAreaCompare';

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
    'rho_pval_win_50ms'
    'rhoRaw_50ms'
    'pvalRaw_50ms'
    'rhoZTrial_50ms'
    'pvalZTrial_50ms'
    'rho_pval_win_150ms'
    'rhoRaw_150ms'
    'pvalRaw_150ms'
    'rhoZTrial_150ms'
    'pvalZTrial_150ms'
    'rho_pval_win_200ms'
    'rhoRaw_200ms'
    'pvalRaw_200ms'
    'rhoZTrial_200ms'
    'pvalZTrial_200ms'
    %'signifRaw_05'
    %'signifRaw_01'
    %'signifZTrial_05'
    %'signifZTrial_01'
    };

t = regexp(useCols,'rhoRaw_(\d*)ms$','tokens');
staticWinSizes = sort(cellfun(@(x) str2double(x{1}),t(~cellfun(@isempty,t))));

%% aggregate required fields into a single table
spkCorrAllTbl = table();
for pa =1:numel(pairAreas)
    pairArea = pairAreas{pa};
    spkCorrAllTbl = [spkCorrAllTbl; spkCorr.(pairArea)(:,useCols)];
end
% recode X_monkey
spkCorrAllTbl.monkey = spkCorrAllTbl.X_monkey;
spkCorrAllTbl.X_monkey = [];


%% Filter data?
filterTxt = {};
outFileFilterSuffix = '';
if filterData
    filters.monk = {'D','E'};
    filt = filters.monk{filterIdx};
    filterTxt = [filterTxt, ['Monk: ' filt]];
    outFileFilterSuffix = [outFileFilterSuffix '_' filt];
    spkCorrAllTbl = spkCorrAllTbl(strcmp(spkCorrAllTbl.monkey,filt),:);
else
    outFileFilterSuffix = '_All';
end

if filterData
    txtFilter = ['Filters : ' char(join(filterTxt',' ; '))]; %#ok<*UNRCH>
else
    txtFilter = '';
end


%% Plot data for each window
for sw = 1:numel(staticWinSizes)
    
    %% For each static window size
    spkCorrAll = spkCorrAllTbl;
    staticWin = staticWinSizes(sw);
    swSuffix = num2str(staticWin,'_%dms');
    outFileSuffix = [outFileFilterSuffix swSuffix]; %#ok<*AGROW>
    rhoRawName = 'rhoRaw';
    spkCorrAll.(rhoRawName) = spkCorrAll.([rhoRawName swSuffix]);
    
    rhoZTrialName = 'rhoZTrial';
    spkCorrAll.(rhoZTrialName) = spkCorrAll.([rhoZTrialName swSuffix]);
    
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
    
    %% get stats
    % Fields for grpstats function:
    grpstatsFields = {'baseCondition','trialType','newPairAreas','condition','alignedName'};
    grpstatsFields = [grpstatsFields';rhoRawName;rhoZTrialName];
    spkCorrAll = spkCorrAll(:,grpstatsFields);
    % separate stats for Positive rho and negative rho
    spkCorrAllStats = struct();
    % Positive Rho
    idx =  spkCorrAll.(rhoRawName)>=0;
    tempStats = grpstats(spkCorrAll(idx,:),...
        {'condition','baseCondition','trialType','alignedName','newPairAreas'},{'mean','std'});
    tempStats.('sem_rhoRaw') = (tempStats.(['std_' rhoRawName]).^2)./sqrt(tempStats.GroupCount);
    tempStats.('sem_rhoZTrial') = (tempStats.(['std_' rhoZTrialName]).^2)./sqrt(tempStats.GroupCount);
    spkCorrAllStats.plusRho = tempStats;
    % Negative Rho
    idx =  spkCorrAll.(rhoRawName)<0;
    tempStats = grpstats(spkCorrAll(idx,:),...
        {'condition','baseCondition','trialType','alignedName','newPairAreas'},{'mean','std'});
    tempStats.('sem_rhoRaw') = (tempStats.(['std_' rhoRawName]).^2)./sqrt(tempStats.GroupCount);
    tempStats.('sem_rhoZTrial') = (tempStats.(['std_' rhoZTrialName]).^2)./sqrt(tempStats.GroupCount);
    spkCorrAllStats.minusRho = tempStats;
    % Absolute Rho
    spkCorrAll.(rhoRawName) = abs(spkCorrAll.(rhoRawName));
    spkCorrAll.(rhoZTrialName) = abs(spkCorrAll.(rhoZTrialName));   
    tempStats = grpstats(spkCorrAll,...
        {'condition','baseCondition','trialType','alignedName','newPairAreas'},{'mean','std'});
    tempStats.('sem_rhoRaw') = (tempStats.(['std_' rhoRawName]).^2)./sqrt(tempStats.GroupCount);
    tempStats.('sem_rhoZTrial') = (tempStats.(['std_' rhoZTrialName]).^2)./sqrt(tempStats.GroupCount);
    spkCorrAllStats.absoluteRho = tempStats;    
    %% Plot it

    pdfFn = [outPdfFn outFileSuffix '.pdf'];
    plotCorrStats(pdfFn,spkCorrAllStats,newPairAreas,txtFilter);
    
end
%% Other functions  

function [] = plotCorrStats(pdfFn,spkCorrAllStats,newPairAreas,txtFilter)
% spkCorrAllStats is a struct
%    .plusRho = all positive rho (including zero)
%    .minusRho = all negative rho
%    .absoluteRho = absolute value of rho

    baseConditions = {'Correct','ErrorChoice','ErrorTiming'};
    alignedNames = {'Baseline','Visual','PostSaccade','PostReward'};
    % find minRho for negative Rho 
    % manually was [-0.2 0.45]
    
    H_plots = getPlotHandles();
    for plusMinus = 1:2
        isPlusRhos = 1; %sum(yMean(~isnan(yMean))>=0)==numel(~isnan(yMean));
        currStats = spkCorrAllStats.plusRho;
        if plusMinus == 2
            currStats = spkCorrAllStats.minusRho;
            isPlusRhos = 0;
        end
        H_plt_Idx = 0;
        for co = 1:numel(alignedNames)
            alignedName = alignedNames{co};
            for ro = 1:numel(baseConditions)
                baseCond = baseConditions{ro};
                tblRowsAccu = strcmp(currStats.alignedName,alignedName) & ...
                    strcmp(currStats.baseCondition,baseCond) & strcmp(currStats.trialType,'Accurate');
                tblRowsFast = strcmp(currStats.alignedName,alignedName) & ...
                    strcmp(currStats.baseCondition,baseCond) & strcmp(currStats.trialType,'Fast');

                % It is possible that some pair-areas or groups may not exist
                % in data.. In these cases, add 'NaN' to data
                %if numel(pairAreasInData)~=numel(newPairAreas)
                tempTbl = table();
                tempTbl.newPairAreas = newPairAreas';
                tempTbl.grpOrder = (1:size(tempTbl,1))';
                accuTbl = outerjoin(tempTbl,currStats(tblRowsAccu,:),'MergeKeys',true);
                fastTbl = outerjoin(tempTbl,currStats(tblRowsFast,:),'MergeKeys',true);
                accuTbl = sortrows(accuTbl,'grpOrder');
                fastTbl = sortrows(fastTbl,'grpOrder');
                % create table for plot 
                pltTbl = table();
                pltTbl.datAreas = newPairAreas';
                pltTbl.datMean = [accuTbl.mean_rhoRaw fastTbl.mean_rhoRaw];
                pltTbl.datStd = [accuTbl.std_rhoRaw fastTbl.std_rhoRaw];
                pltTbl.datSem = [accuTbl.sem_rhoRaw fastTbl.sem_rhoRaw];
                pltTbl.datCounts = [accuTbl.GroupCount fastTbl.GroupCount];
                 
                % plot it
                H_plt_Idx = H_plt_Idx + 1;
                plotGroupedBarsWithErrBars(H_plots(H_plt_Idx),pltTbl.datAreas,...
                    pltTbl.datMean,pltTbl.datStd,pltTbl.datSem,pltTbl.datCounts,isPlusRhos);
            end
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


function [] = plotGroupedBarsWithErrBars(H_axes,xData,yMean,yStd,ySem,nPairs,isPlusRhos) %#ok<*INUSL>

    % find minRho for negative Rho 
    % manually was [-0.2 0.45]
    yLims = [-0.4 0.4];
    useSem = 1;
    axes(H_axes);
    accColor = [1 0 0];
    fastColor = [0 0.5 0];
    
    %x = categorical(xData(:,1));
    x = 1:numel(xData(:,1));
    hBar = bar(x,yMean,'FaceAlpha',0.4,'BarWidth',0.90);
    set(hBar(1),'FaceColor',accColor); set(hBar(2),'FaceColor',fastColor);
    set(gca,'xticklabels',xData(:,1),'TickLabelInterpreter','none','XTickLabelRotation',20);
    set(gca,'YGrid','on')
    ylabel('Mean r_{sc} \pm SEM (abs.)');
    set(gca,'Ylim',yLims);
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
    he = errorbar(ctr,yMean',yErr','Marker','o','MarkerSize',5,'LineStyle',...
        'none','LineWidth',0.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
    % add text values to plot
    if isPlusRhos
        set(he,'YPositiveDelta',[]);
        y = 0.2;
        va = 'bottom';
    else
        set(he,'YNegativeDelta',[]);
         y = -0.2;
        va = 'top';
   end    
    
    
    [txtPairArea,txtAcc,txtFast] = arrayfun(@(x) deal(sprintf('%14s',xData{x,1}),...
        sprintf('%0.2f\\pm%0.2f %4d',yMean(x,1),yErr(x,1),nPairs(x,1)),...
        sprintf('%0.2f\\pm%0.2f %4d',yMean(x,2),yErr(x,2),nPairs(x,2))),...
        (1:size(yMean,1))','UniformOutput',false);
    txtPairArea = [sprintf('%14s','AREA');txtPairArea];
    txtAcc = [sprintf('%14s','\mu\pmSEM   n');txtAcc];
    txtFast = [sprintf('%14s','\mu\pmSEM   n');txtFast];
%     


    h_pa = text(ctr(2,2)+0.5,y,txtPairArea,'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',7);
    ext = get(h_pa,'Extent');
    pos = get(h_pa,'Position');
    h_txt1 = text(pos(1)+ext(3),y,txtAcc,'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',7,'color',accColor);
    ext = get(h_txt1,'Extent');
    pos = get(h_txt1,'Position');
    h_txt2 = text(pos(1)+ext(3)+0.1,y,txtFast,'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',7,'color',fastColor);
    
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
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount));
        end
    end


end
    
    

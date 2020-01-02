% load extracted data
dat = load('dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticNew.mat');

outDir = 'dataProcessed/analysis/spkCorr/summary';
fns = fieldnames(dat);
% area pairs to plot (in-order). Do not change the order below
pairAreas = {
    'SEF_SEF'
    'FEF_FEF'
    %'SC_SC'
    };
assert(numel(pairAreas)== 2|numel(pairAreas)== 3 ,...
    sprintf('The number of area-pairs [%d] must be [2 or 3]',numel(pairAreas)));
colNames = dat.(fns{1}).Properties.VariableNames';
t = regexp(colNames,'rhoRaw_(\d*)ms$','tokens');
staticWinSizes = sort(cellfun(@(x) str2double(x{1}),t(~cellfun(@isempty,t))));

plotUnsignedCorr = 0;
saveFigFlag = 1;

if plotUnsignedCorr
    absOrSigned = '(abs.)';
    fnSuffix = '_Absolute';
else
    absOrSigned = '';
    fnSuffix = '';
end

for sw = 1:numel(staticWinSizes)
    
    %% For each static window size
    staticWin = staticWinSizes(sw);
    swSuffix = num2str(staticWin,'_%dms');
    rhoPvalWinName = ['rho_pval_win' swSuffix];
    rhoName = ['rhoRaw' swSuffix];
    
    [alignedNames,idx] = unique(dat.(pairAreas{1}).alignedName,'stable');
    conditions = unique(dat.(pairAreas{1}).condition);
    alignedOn = dat.(pairAreas{1}){idx,{'alignedEvent'}};
    rscTimeWins =  dat.(pairAreas{1}){idx,{rhoPvalWinName}};
    alignedOnTimeWinsStr = cellfun(@(x,y) [sprintf('Aligned on: %s',x) ...
        ', Spike Counts Win:[' num2str(y,'%d ') '], Width:[' swSuffix(2:end) '] '],...
        alignedOn,rscTimeWins,'UniformOutput',false);
    % there will be nAlignedNames figures - Each figure will have
    % SEF-SEF,FEF-FEF, SC-SC pair correlations plotted againts distance between pairs
    % for all conditions: (Fast / Accurate (Correct, ErrorChoice, ErrorTiming))
    % for epochs (fig1) Baseline, (fig2) Visual, (fig3) PostSaccade, (fig4)
    % PostReward
    warning('off')
    
    for figNo = 1:numel(alignedOn)
        epoch = alignedNames{figNo};
        figTitleStr = [epoch ' : ' alignedOnTimeWinsStr{figNo} absOrSigned];
        fn = fullfile(outDir,['corrByDistDistrib_' epoch fnSuffix swSuffix '.pdf']);
        [aspectRatio,H_figure] = getFigHandle();
        [H_plots] = getPlotHandlesB(H_figure, numel(pairAreas));
        annotation('textbox','Position',[0.05 0.965 0.80 0.05],'String',figTitleStr,...
            'Interpreter','none','EdgeColor','none','FontWeight','bold',...
            'FitBoxToText','on','HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontSize',16);
        H_plot_Idx = 0;
        % col1: {AccurateCorrect-SEF,FEF,SC; FastCorrect-SEF,FEF,SC
        plotColConds = {'Correct','ErrorChoice','ErrorTiming'};
        titleColors ={'g','r'};% Fast/Accurate
        for co = 1:numel(plotColConds)
            plotRowConds = {'Fast','Accurate'};
            plotColName = plotColConds{co};
            for ro = 1:numel(plotRowConds)
                titleColor = titleColors{ro};
                plotRowName = plotRowConds{ro};
                plotPairAreas = strrep(pairAreas,'_','-');
                for ar = 1:numel(plotPairAreas)
                    plotPairArea = plotPairAreas{ar};
                    currDat = dat.(strrep(plotPairArea,'-','_'));
                    filteredFlag = strcmp(currDat.alignedName,epoch) ...
                        & strcmp(currDat.condition,[plotRowName plotColName]) ...
                        & strcmp(currDat.pairAreas,plotPairArea);
                    currDat = currDat(filteredFlag,:);
                    currDat = currDat(~isnan(currDat.XY_Dist),:);
                    
                    if plotUnsignedCorr
                        currDat.rhoRaw = abs(currDat.(rhoName));
                    else
                        currDat.rhoRaw = currDat.(rhoName);
                    end
                    
                    % (1) percentage of pairs vs binned-spikecorr,
                    % (2) spike corr vs binned-distance
                    H_plot_Idx = H_plot_Idx + 1;
                                        
                    plotRhoDistrib(H_plots(H_plot_Idx),currDat.rhoRaw);                                       

                    ylabel(['Percentage of ' plotPairArea ' Pairs'],'FontWeight','bold','FontAngle','italic');
                    xlabel('Spike Count Correlation (r_{sc})','Interpreter','tex','FontWeight','bold','FontAngle','italic');
                    titlePos = get(gca,'Position');
                    
                    % (2) spike corr vs binned-distance
                    H_plot_Idx = H_plot_Idx + 1;
                    plotRhoVsDistDistrib(H_plots(H_plot_Idx),currDat.XY_Dist,currDat.rhoRaw);
                    ylabel(['\rho ' absOrSigned plotPairArea ' pairs'],'Interpreter','tex','FontWeight','bold','FontAngle','italic');
                    xlabel('Distance between units (mm)','Interpreter','tex','FontWeight','bold','FontAngle','italic');
                    
                    if ar == 1
                        thPos(1) = titlePos(1) + 0.01 ;
                        thPos(2) = titlePos(2) + titlePos(4);
                        thPos(3) = 0.2;
                        thPos(4) = 0.03;
                        th = annotation('textbox','String',[plotRowName plotColName],'Position', thPos,...
                            'Color',titleColor,'FontSize',14,'FontWeight','bold','FontAngle','italic',...
                            'HorizontalAlignment','center','EdgeColor','none');
                    end
                end
            end
            drawnow;
        end
        if saveFigFlag
            saveFigPdf(fn);
            delete(gcf)
        end
    end
end

function [] = plotRhoVsDistDistrib(H_axes,rawxDat,rawyDat)
    rawxDat = double(rawxDat);
    rawyDat = double(rawyDat);
    idx = ~isnan(rawxDat);
    distBins = [0.5:1:8];
    xLims = [-1 8];
    xTicks = [0:7];
    yLims = [-0.3 0.3];
  
    datTest = table();
    datTest.rho = rawyDat(idx);
    datTest.dist = rawxDat(idx);
    datTest.newDist = rawxDat(idx);

    for r = 1:numel(distBins)-1
        lo = distBins(r);
        hi = distBins(r+1);
        idx = datTest.dist>lo & datTest.dist<=hi;
        datTest.newDist(idx) = (lo + hi)/2;
    end

    plusIdx = datTest.rho>=0;
    minusIdx = datTest.rho<0;
    datStats = struct();
  
    plusRho = grpstats(datTest(plusIdx,:),{'newDist'},{'mean','std'});
    plusRho.sem_rho = (plusRho.std_rho.^2)./sqrt(plusRho.GroupCount);
  
    minusRho = grpstats(datTest(minusIdx,:),{'newDist'},{'mean','std'});
    minusRho.sem_rho = (minusRho.std_rho.^2)./sqrt(minusRho.GroupCount);
       
     axes(H_axes)
     % do plus rho vs distance
     x = plusRho.newDist;
     y = plusRho.mean_rho;
     yerr = plusRho.sem_rho;
     ycounts = plusRho.GroupCount;
     h1 = bar(x,y,'BarWidth',1,'FaceAlpha',0.5);
     
     set(gca,'XLim',xLims,'XTick',xTicks);
     set(gca,'YLim',yLims);
     
     grid on
     set(gca,'YMinorGrid','on');
     
     hold on
     he1 = errorbar(x,y,yerr,'YPositiveDelta',[],'Marker','o','MarkerSize',6,'LineStyle',...
         'none','LineWidth',1.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','c');
     %set(gca,'YLim',[0 max(y)*1.1])
     muStdN = num2str(ycounts,'%d');
     htxt1 = text(x+0.1,y+0.005,muStdN,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8);

     % do minus rho vs distance
     x = minusRho.newDist;
     y = minusRho.mean_rho;
     yerr = minusRho.sem_rho;
     ycounts = minusRho.GroupCount;
     h2 = bar(x,y,'BarWidth',1,'FaceAlpha',0.5);
     he2 = errorbar(x,y,yerr,'YNegativeDelta',[],'Marker','o','MarkerSize',6,'LineStyle',...
         'none','LineWidth',1.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','c');
     %set(gca,'YLim',[0 max(y)*1.1])
     muStdN = num2str(ycounts,'%d');
     htxt2 = text(x+0.1,y-0.005,muStdN,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
     
end

function [] = plotRhoDistrib(H_axes,rawDat)
    
     rawDat = double(rawDat);
     minMaxRho = round(minmax(rawDat(:)'),1); % to 1st decimal
     maxRho = max(abs(minMaxRho));
     maxRho = 0.4;
     binWidth = 0.02;
     plusBins = [0:binWidth:maxRho];
     minusBins = [0:-binWidth:-maxRho];
     recodeRho = rawDat;
     for r = 1:numel(plusBins)-1
         lo = plusBins(r);
         hi = plusBins(r+1);
         idx = recodeRho>lo & recodeRho<=hi;
         recodeRho(idx) = (lo + hi)/2;
     end
     for r = 1:numel(minusBins)-1
         hi = minusBins(r);
         lo = minusBins(r+1);
         idx = recodeRho>=lo & recodeRho<hi;
         recodeRho(idx) = (lo + hi)/2;
     end
     
     bEdge = [-maxRho:binWidth:maxRho];
     y = histcounts(recodeRho,bEdge);
     y(y==0) = NaN;
     %x = bEdge(1:end-1)+binWidth/2;     
     x = bEdge(2:end);
     xMinor = [min(x):binWidth:max(x)];
     
     xPlus = x>0;
     xMinus = x<0;
     yPercent = y.*100/nansum(y);
     
     maxy = max(yPercent)*1.2;
     yLims = [-maxy maxy];
     
     axes(H_axes)
     % plus plot
     hplus = bar(x(xPlus),yPercent(xPlus),'BarWidth',1,'FaceAlpha',0.5);
     hold on
     h1 = text(x(xPlus),yPercent(xPlus),num2str(y(xPlus)','%d'),...
         'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8); 
     ylim(yLims)
     grid on
     set(gca,'XMinorGrid','on')
     yrange = range(get(gca,'ylim'));
     % annotate
     idx = rawDat>0;
     scatter(nanmean(rawDat(idx)),maxy*0.9,100,'v','filled','MarkerEdgeColor','r','MarkerFaceColor','r'); 
     line([nanmean(rawDat(idx)),nanmean(rawDat(idx))],[0 maxy*0.9],'color','r','LineStyle','--','LineWidth',1);
     text(nanmean(rawDat(idx))+0.02,maxy*0.9,sprintf('%0.2f (%d)',nanmean(rawDat(idx)),numel(rawDat(idx))),'FontSize',10,'FontWeight','bold');
     % minus plot
     hminus = bar(-x(xMinus),-yPercent(xMinus),'BarWidth',1,'FaceAlpha',0.5);
     h2 = text(-x(xMinus),-yPercent(xMinus),num2str(y(xMinus)','%d'),...
         'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8); 
     idx = rawDat<0;
     scatter(-nanmean(rawDat(idx)),-maxy*0.9,100,'^','filled','MarkerEdgeColor','r','MarkerFaceColor','r'); 
     line([-nanmean(rawDat(idx)),-nanmean(rawDat(idx))],[0 -maxy*0.9],'color','r','LineStyle','--','LineWidth',1);
     text(-nanmean(rawDat(idx))+0.02,-maxy*0.9,sprintf('%0.2f (%d)',nanmean(rawDat(idx)),numel(rawDat(idx))),'FontSize',10,'FontWeight','bold');

end
    
%%
function [aspectRatio, H_Figure] = getFigHandle()
    H_Figure = newFigure();
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);
end

function  [H_plots] = getPlotHandlesB(H_Figure,nAreas)
    % number of plot rows per column = nAreas for fast, nAreas for accurate
    % Hence nRows = nAreas * 2 (example, plot SEF-SEF,FEF-FEF,SC-SC)
    % A 4 by 3 grid; Each Cell will have 1 plot
    % Or a 6 by 3 grid; Each Cell will have 1 plot

    nRows = nAreas*2;
    if nRows == 6
        pltH = 0.11;
        pltW = 0.26;
        gutter = 0.03;
        offsetsX = 0.05:(pltW+gutter*2):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 6 row-starts
        offsetsY(4:end) = offsetsY(4:end) - gutter;
    elseif nRows == 4
        pltH = 0.17;
        pltW = 0.275;
        gutter = 0.05;
        offsetsX = 0.05:(pltW+gutter):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 4 row-starts
        offsetsY(3:end) = offsetsY(3:end) - gutter;
    end

    pltCount = 0;
    for cols = 1:3
        pltW2 = (pltW-gutter)/2;
        for ros = 1:nRows
            for subp = 1:2
                % (1) percentage of pairs vs binned-spikecorr,
                % (2) spike corr vs binned-distance
                pos(1) = offsetsX(cols) + (pltW2+gutter/1.5)*(subp-1);
                pos(3:4) = [pltW2 pltH];
                pos(2) = offsetsY(ros);
                pltCount = pltCount + 1;
                H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
            end
        end
    end

end


function saveFigPdf(fn)
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',screenposition(3:4),...
        'PaperOrientation','landscape');
    fprintf('Saving figure to: %s\n',fn);
    print(fn,'-dpdf','-opengl')
    drawnow
end


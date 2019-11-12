% load extracted data
%dat = load('dataProcessed/analysis/spkCorr/spkCorrAllPairsStatic.mat');
plotUnsignedCorr = 1;
saveFigFlag = 1;
outDir = 'dataProcessed/analysis/spkCorr';
fns = fieldnames(dat);
% area pairs to plot (in-order). Do not change the order below
pairAreas = {
    'SEF_SEF'
    'FEF_FEF'
    %'SC_SC'
    };
assert(numel(pairAreas)== 2|numel(pairAreas)== 3 ,...
    sprintf('The number of area-pairs [%d] must be [2 or 3]',numel(pairAreas))); 

[alignedNames,idx] = unique(dat.(pairAreas{1}).alignedName,'stable');
conditions = unique(dat.(pairAreas{1}).condition);
alignedOn = dat.(pairAreas{1}){idx,{'alignedEvent'}};
rscTimeWins =  dat.(pairAreas{1}){idx,{'rho_pval_win'}};
alignedOnTimeWinsStr = cellfun(@(x,y) ['[' sprintf('Aligned on: %s',x) ', Spike Counts Win: [' num2str(y,'%d  ') ']]'],...
                       alignedOn,rscTimeWins,'UniformOutput',false);
% there will be nAlignedNames figures - Each figure will have 
% SEF-SEF,FEF-FEF, SC-SC pair correlations plotted againts distance between pairs
% for all conditions: (Fast / Accurate (Correct, ErrorChoice, ErrorTiming))
% for epochs (fig1) Baseline, (fig2) Visual, (fig3) PostSaccade, (fig4)
% PostReward
warning('off')
suff = ' ';
if plotUnsignedCorr
    suff = '(abs.) ';
end
for figNo = 1:numel(alignedOn)
    epoch = alignedNames{figNo};
    figTitleStr = [epoch ' : ' alignedOnTimeWinsStr{figNo} suff];
    [aspectRatio,H_figure] = getFigHandle();
    [H_plots] = getPlotHandlesB(H_figure, numel(pairAreas));
    annotation('textbox','Position',[0.05 0.97 0.80 0.05],'String',figTitleStr,...
               'Interpreter','tex','EdgeColor','none','FontWeight','bold',...
               'FitBoxToText','on','HorizontalAlignment','center',...
               'VerticalAlignment','baseline','FontSize',18);    
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
                  currDat.rhoRaw = abs(currDat.rhoRaw);
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
                  ylabel(['\rho ' suff plotPairArea ' pairs'],'Interpreter','tex','FontWeight','bold','FontAngle','italic');
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
         suff2 = '';
         if plotUnsignedCorr
             suff2 = '_Unsigned';
         end
         fn = fullfile(outDir,['Summary_Distrib_' epoch '_SpkCorr' suff2 '.pdf']);
         saveFigPdf(fn);
         delete(gcf)
     end
end

function [] = plotRhoVsDistDistrib(H_axes,rawxDat,rawyDat)
     rawxDat = double(rawxDat);
     rawyDat = double(rawyDat);
     minMaxDist = round(minmax(rawxDat(:)'),1); % to 1st decimal
     binWidth = 1.0;%mm
     %bins = (-binWidth:binWidth:minMaxDist(2))'; 
     bins = (minMaxDist(1):binWidth:ceil(minMaxDist(2)))';
     y = arrayfun(@(x) mean(rawyDat(rawxDat>=bins(x) & rawxDat<bins(x+1))), (1:numel(bins)-1)');
     ystd = arrayfun(@(x) std(rawyDat(rawxDat>=bins(x) & rawxDat<bins(x+1))), (1:numel(bins)-1)');
     ycounts = arrayfun(@(x) sum(rawxDat>=bins(x) & rawxDat<bins(x+1)), (1:numel(bins)-1)');
     x = bins(1:end-1);
     
     axes(H_axes)
     h = bar(x,y,'BarWidth',1,'FaceAlpha',0.5);
     hold on
     he = errorbar(x,y,ystd,'Marker','o','MarkerSize',10,'LineStyle',...
         'none','LineWidth',1.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','c')
     %set(he,'YPositiveDelta',[])
     set(gca,'YLim',[0 max(y)*1.1])
     mustdn = [num2str(ycounts(:),'%d')]
     h1 = text(x+0.1,y+0.005,mustdn,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8);
     
end

function [] = plotRhoDistrib(H_axes,rawDat)
    
     rawDat = double(rawDat);
     minMaxRho = round(minmax(rawDat(:)'),1); % to 1st decimal
     binWidth = 0.05;
     bins = minMaxRho(1):binWidth:minMaxRho(2);
     y = histcounts(rawDat,bins);
     y(y==0) = NaN;
     ypercent = y.*100/nansum(y);
     x = bins(1:end-1) + binWidth/2;
     
     axes(H_axes)
     h = bar(x,ypercent,'BarWidth',1,'FaceAlpha',0.5);
     hold on
     h1 = text(double(x),ypercent,num2str(y(:),'%d'),...
         'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8); 
     ylim([0 max(ypercent)+5])
     yrange = range(get(gca,'ylim'));
     % annotate
     scatter(nanmean(rawDat),yrange*0.9,100,'v','filled','MarkerEdgeColor','r','MarkerFaceColor','r'); 
     line([nanmean(rawDat),nanmean(rawDat)],[0 yrange*0.9],'color','r','LineStyle','--','LineWidth',1);
     text(nanmean(rawDat),yrange*0.95,sprintf('%0.2f (%d)',nanmean(rawDat),numel(rawDat)),'FontSize',12,'FontWeight','bold');
end
    
%%
function [aspectRatio, H_Figure] = getFigHandle()
%set(0, 'DefaultFigureColormap', jet(64));
set(0,'units','pixels');
set(0,'defaulttextfontsize',6,...
    'defaulttextfontname','Arial',...
    'defaultaxesfontsize',6,...
    'defaultaxeslinewidth',0.05);
margin = 10; %pixels
%ss=get(0,'ScreenSize');
% optimized for this size on my macbookpro
ss = [1 1 1680 1050];
ss = [20 20 1500 990];
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
end

function  [H_plots] = getPlotHandlesB(H_Figure,nAreas)
    % number of plot rows per column = nAreas for fast, nAreas for accurate
    % Hence nRows = nAreas * 2 (example, plot SEF-SEF,FEF-FEF,SC-SC)
    % A 4 by 3 grid; Each Cell will have 1 plot
    % Or a 6 by 3 grid; Each Cell will have 1 plot

    nRows = nAreas*2;
    if nRows == 6
        pltH = 0.12;
        pltW = 0.26;
        gutter = 0.03;
        offsetsX = 0.05:(pltW+gutter*2):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 6 row-starts
        offsetsY(4:end) = offsetsY(4:end) - gutter;
    elseif nRows == 4
        pltH = 0.18;
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


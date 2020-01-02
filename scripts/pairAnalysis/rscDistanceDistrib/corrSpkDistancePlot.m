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
    signif05Name = ['signifRaw_05' swSuffix];
    signif01Name = ['signifRaw_01' swSuffix];
    
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
        fn = fullfile(outDir,['corrByDistance_' epoch fnSuffix swSuffix '.pdf']);
        [aspectRatio,H_figure] = getFigHandle();
        [H_plots] = getPlotHandles(H_figure, numel(pairAreas));
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
                    
                    % round XY_Dist to nearest microns
                    roundToMs = 200;
                    currDat.XY_DistBinned = round(currDat.XY_Dist*1000/roundToMs).*(roundToMs/1000);
                    % 
                    plusDat = currDat(currDat.rhoRaw>=0,:);
                    minusDat = currDat(currDat.rhoRaw<0,:);                    
                                        
                    plusStats = grpstats(plusDat(:,{'XY_DistBinned','rhoRaw'}),'XY_DistBinned',{'mean','std'});
                    minusStats = grpstats(minusDat(:,{'XY_DistBinned','rhoRaw'}),'XY_DistBinned',{'mean','std'});
                   
                    % plot spike corr scatter
                    H_plot_Idx = H_plot_Idx + 1;
                    % Positive correlations
                    [legTxt1] = plotRscScatterStats(H_plots(H_plot_Idx),plusDat.XY_DistBinned,plusDat.rhoRaw,...
                        plusStats.XY_DistBinned,plusStats.mean_rhoRaw,plusDat.(signif05Name),plusDat.(signif01Name));
                    % Negative correlations
                    [legTxt2] = plotRscScatterStats(H_plots(H_plot_Idx),minusDat.XY_DistBinned,minusDat.rhoRaw,...
                        minusStats.XY_DistBinned,minusStats.mean_rhoRaw,minusDat.(signif05Name),minusDat.(signif01Name));
                    % add legend for positive and negative corr
                    leg1 = regexprep(legTxt1,' pairs.*','(+) ');
                    leg2 = regexprep(legTxt2,' pair','(-) pair');
                    legend(strcat(leg1,'  ',leg2),'Location','northeast','Interpreter','tex','Box','off');
                    
                    ylabel(['\rho ' absOrSigned plotPairArea ' pairs'],'Interpreter','tex','FontWeight','bold','FontAngle','italic');
                    xlabel('Distance between units (mm)','Interpreter','tex','FontWeight','bold','FontAngle','italic');
                    
                    if ar==1
                        title([plotRowName plotColName],'Color',titleColor,'FontSize',14,'FontWeight','bold','FontAngle','italic');
                    end
                end
            end
            drawnow;
        end
        % scale Ylims of all plots to be same
            
        if saveFigFlag
            saveFigPdf(fn);
            delete(gcf)
        end
    end
end

function [legTxt] = plotRscScatterStats(H_axes,x,y,xmu,ymu,sig05,sig01)
% x = x values for y = raw spkcorr-rho value
% xmu = x values for ymu = mean spkcorr-rho value
% ystd = std of raw spkcorr-rho values (for plotting +/-std around mean
% sig05 = boolean, raw spkcorr values <= pval of 0.05
% sig01 = boolean, raw spkcorr values <= pval of 0.01
%    yScaleFixed = false;
%    yLims = [-1 1]; % set to min/max for Pearson corr-coeff
%     if ~yScaleFixed
%         yLims = minmax(y(:)');% row
%         % round to 2 decimals
%         yLims = round((yLims + [-0.05 +0.05]),2);
%     end

    if isempty(x) || isempty(y)
        return;
    end
    
    col_k = [0 0 0]; col_m = [1 0 1]; col_c = [0 1 1]; col_b=[0 0 1]; col_r = [1 0 0];
    axes(H_axes);
    notSignif = ~or(sig01,sig05);
    scatter(x(notSignif),y(notSignif),'o','filled','MarkerFaceColor',col_k,'MarkerFaceAlpha',0.4)
    hold on
    xlim([min(x)-0.5 max(x)+0.5])
    scatter(x(sig05),y(sig05),'o','filled','MarkerFaceColor',col_b,'MarkerFaceAlpha',0.4)
    scatter(x(sig01),y(sig01),'o','filled','MarkerFaceColor',col_r,'MarkerFaceAlpha',0.4)
    % fit linear model
    mdl = fitlm(x,y,'linear');
    [ypred,ci05] = predict(mdl,xmu(:),'Alpha',0.05);
    plot(xmu,ypred,'Color',col_k,'LineWidth',1.5,'HandleVisibility','off');
    fill([xmu(:)' fliplr(xmu(:)')],[ci05(:,1)' fliplr(ci05(:,2)')],col_k,'FaceAlpha',0.1,'LineStyle','none','HandleVisibility','off');
    %plot(xmu(:),ci05,'Color',col_b,'HandleVisibility','off');
    % draw the mean scatter
    scatter(xmu,ymu,60,'Marker','diamond','LineWidth',1,'MarkerEdgeColor',col_k)
    plot(xmu,ymu,'Marker','none','Color',col_k,'LineStyle','--','HandleVisibility','off');
    % Draw zero line
    line(get(gca,'XLim'),[0 0]);
    
    % annotate counts for each plot
    legTxt = {sprintf('%d pairs p>0.05',sum(notSignif)),...
        sprintf('%d pairs p<=0.05',sum(sig05)),...
        sprintf('%d pairs p<=0.01',sum(sig01)),...
        sprintf('%d pairs mean corr.',numel(y))...
        };
    % legend(legTxt,'Location','northeast','Interpreter','tex','Box','off');
     
end
    
    
%%
function [aspectRatio, H_Figure] = getFigHandle()
    H_Figure = newFigure();
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);
end

function  [H_plots] = getPlotHandles(H_Figure,nAreas)
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
        gutter = 0.04;
        offsetsX = 0.05:(pltW+gutter):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 4 row-starts
        offsetsY(3:end) = offsetsY(3:end) - gutter;
    end   
    
    pltCount = 0;
    for cols = 1:3
        pos(1) = offsetsX(cols);
        pos(3:4) = [pltW pltH];
        for ros = 1:nRows
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
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


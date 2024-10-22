% load extracted data
dat = load('dataProcessed/analysis/spkCorr_TimeWin_Z/rscSignalNoiseStatic.mat');
saveFigFlag = 1;
outDir = 'dataProcessed/analysis/spkCorr_TimeWin_Z';
fns = fieldnames(dat);
% area pairs to plot (in-order). Do not change the order below
pairAreas = {
    'SEF_SEF'
    'SEF_FEF'
    'SEF_SC'
    'SEF_NSEFN'
    'FEF_FEF'
    'FEF_SC'
    'FEF_NSEFN'
    'SC_SC'
    'SC_NSEFN'
    'NSEFN_NSEFN' 
    };
conditions = unique(dat.(pairAreas{1}).condition);
[alignedNames,idx] = unique(dat.(pairAreas{1}).alignedName,'stable');
alignedOn = dat.(pairAreas{1}){idx,{'alignedEvent'}};
rscTimeWins =  dat.(pairAreas{1}){idx,{'rho_pval_win_Z'}};
alignedOnTimeWinsStr = cellfun(@(x,y) ['[' sprintf('Aligned on: %s',x) ', TimeWin: [' num2str(y,'%d  ') ']]'],...
                       alignedOn,rscTimeWins,'UniformOutput',false);
% there will be nConditions * nAlignedNames figures
warning('off')
for cond = 1:numel(conditions)
    condition = conditions{cond};  
    for an = 1:numel(alignedNames)
        epoch = alignedNames{an};
        [aspectRatio,H_figure] = getFigHandle();
        [H_plots] = getPlotHandles(H_figure);
        H_plot_Idx = 0;
        for pa = 1:numel(pairAreas)
            pairArea = pairAreas{pa};
            currDat = dat.(pairArea);
            currDat = currDat(~isnan(currDat.XY_Dist),:);
            filteredFlag = strcmp(currDat.condition,condition) ...
                & strcmp(currDat.alignedName,epoch)...
                & strcmp(currDat.pairAreas,strrep(pairArea,'_','-'));
            currDat = currDat(filteredFlag,:);
            % round XY_Dist to nearest 100 microns
            roundToMs = 200;
            currDat.XY_DistBinned = round(currDat.XY_Dist*1000/roundToMs).*(roundToMs/1000);
            currDatStats = grpstats(currDat(:,{'XY_DistBinned','rhoSignal','rhoNoise'}),'XY_DistBinned',{'mean','std'});

            % plot signal corr
            H_plot_Idx = H_plot_Idx + 1;            
            plotRscScatterStats(H_plots(H_plot_Idx),currDat.XY_DistBinned,currDat.rhoSignal,...
                                currDatStats.XY_DistBinned,currDatStats.mean_rhoSignal,currDatStats.std_rhoSignal,...
                                currDat.signifSignal_05,currDat.signifSignal_01);
            ylabel('\rho Signal','Interpreter','tex'); 
            % plot noise corr
            H_plot_Idx = H_plot_Idx + 1;            
            plotRscScatterStats(H_plots(H_plot_Idx),currDat.XY_DistBinned,currDat.rhoNoise,...
                                currDatStats.XY_DistBinned,currDatStats.mean_rhoNoise,currDatStats.std_rhoNoise,...
                                currDat.signifNoise_05,currDat.signifNoise_01);
            ylabel('\rho Noise','Interpreter','tex'); 
            xlabel('Distance between units (mm)')
            %title(sprintf('%s, %s, %s',pairArea,condition,epoch),'Interpreter','none')
        end
        a = annotation('textbox',[0.25 0.95 0.5 0.05],'String',[condition '  -  ' epoch '   ' alignedOnTimeWinsStr{an}],...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontSize',30,'FontWeight','bold','FontAngle','italic','FitBoxToText','on',...
            'EdgeColor','none','Interpreter','none');
        if saveFigFlag
        fn = fullfile(outDir,['Summary_' condition '_' epoch '.pdf']);
        saveFigPdf(fn);
        delete(gcf)
        end
    end
    
end


function [] = plotRscScatterStats(H_axes,x,y,xmu,ymu,ystd,sig05,sig01)
    if isempty(x) | isempty(y)
        return;
    end
    
    col_k = [0 0 0]; col_m = [1 0 1]; col_c = [0 1 1]; col_b=[0 0 1]; col_r = [1 0 0];
    axes(H_axes);

    scatter(x,y,'o','filled','MarkerFaceColor',col_k,'MarkerFaceAlpha',0.3)
    hold on
    xlim([min(x)-0.5 max(x)+0.5])
    ylim([-1 1])
    scatter(x(sig05),y(sig05),'o','filled','MarkerFaceColor',col_b,'MarkerFaceAlpha',0.4)
    scatter(x(sig01),y(sig01),'o','filled','MarkerFaceColor',col_r,'MarkerFaceAlpha',0.4)
    % fit linear model
    mdl = fitlm(x,y,'linear');
    [ypred,ci05] = predict(mdl,xmu(:),'Alpha',0.05);
    plot(xmu,ypred,'Color',col_k,'LineWidth',1.5);
    fill([xmu(:)' fliplr(xmu(:)')],[ci05(:,1)' fliplr(ci05(:,2)')],col_k,'FaceAlpha',0.1,'LineStyle','none','HandleVisibility','off');
    %plot(xmu(:),ci05,'Color',col_b,'HandleVisibility','off');
    % draw the mean scatter
    plot(xmu,ymu,'d','Color',col_k,'LineStyle','--','HandleVisibility','off');
    
    hold off
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

function  [H_plots] = getPlotHandles(H_Figure)
    % a 4 by 4 grid; Only the tril is populated
    % each cell will have 2 plots    
    offsetsX =[0.05 (1:3).*0.23+0.05]; % for 4 column-starts
    offsetsY = 0.85:-0.23:0.05; % for 4 row-starts  
    pltW = 0.20;
    pltH = 0.08;
    gutter = 0.03;
    pltCount = 0;
    areaLabels = {'SEF','FEF','SC','NSEFN'};
    for cols = 1:4
        pos(1) = offsetsX(cols);        
        pos(3:4) = [pltW pltH];
        % add area label for cols
        coLabs(cols) = annotation('textbox',[offsetsX(cols)+0.05 offsetsY(cols)+0.07 0.10 0.05],'String',areaLabels{cols},...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontSize',18,'FontWeight','bold','FontAngle','italic','FitBoxToText','on',...
            'EdgeColor','none','Interpreter','none');

        for ros = cols:4
            % signal plot
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
            % noise plot
            pos(2) = offsetsY(ros) - pltH - gutter;
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
            if cols == 1
                % add area labels for rows
                roLabs(ros) = annotation('textbox',[0.005 pos(2)+pos(4)-0.01 0.10 0.05],'String',areaLabels{ros},...
                    'HorizontalAlignment','left','VerticalAlignment','middle',...
                    'FontSize',18,'FontWeight','bold','FontAngle','italic','FitBoxToText','on',...
                    'EdgeColor','none','Interpreter','none');
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
    print(fn,'-dpdf','-painters')
    drawnow
end


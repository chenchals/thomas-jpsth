% load extracted data
%dat = load('dataProcessed/analysis/spkCorr/rscSignalNoiseStatic.mat');
fns = fieldnames(dat);
% area pairs to plot (in-order)
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
rscTimeWins =  dat.(pairAreas{1}){idx,{'rho_pval_win_Z'}};

for cond = 1:numel(conditions)
    condition = conditions{cond};
    for an = 1:numel(alignedNames)
        epoch = alignedNames{an}  ;    
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
            % plot signal
            x = currDat.XY_DistBinned;
            y = currDat.rhoSignal;
            sig05 = currDat.signifSignal_05;
            sig01 = currDat.signifSignal_01;
            
            xmu = currDatStats.XY_DistBinned';
            ymu = currDatStats.mean_rhoSignal';
            ystd = currDatStats.std_rhoSignal';
            plotRscScatterStats(x,y,xmu,ymu,ystd,sig05,sig01);
            title(sprintf('%s, %s, %s',pairArea,condition,epoch),'Interpreter','none')
        end
    end
    
end


function [] = plotRscScatterStats(x,y,xmu,ymu,ystd,sig05,sig01)
    gca;
    col_k = [0 0 0]; col_m = [1 0 1]; col_c = [0 1 1]; col_b=[0 0 1]; col_r = [1 0 0];
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
    fill([xmu fliplr(xmu)],[ci05(:,1)' fliplr(ci05(:,2)')],col_k,'FaceAlpha',0.1,'LineStyle','none','HandleVisibility','off');
    %plot(xmu(:),ci05,'Color',col_b,'HandleVisibility','off');
    % draw the mean scatter
    plot(xmu,ymu,'d','Color',col_k,'LineStyle','--','HandleVisibility','off');
    
    hold off
end
    
    
%%
function [H_Figure] = getFigHandle()
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
end

function saveFigAs(fn)

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


function [] = corrSpkCountPlot(spkCountFile,pdfOutputDir,savePdfFlag)
%CORRSPKCOUNTPLOT Summary of this function goes here

%% Put this in the function....
smoothBinWidthMs = 5;
fx_vecSmooth = @(x,w) smoothdata(x,'movmean',w,'omitnan');

conditionPairs = {
    {'FastErrorChoice','AccurateErrorChoice'};
    {'FastErrorTiming','AccurateErrorTiming'};
    };
pdfPrefixMap = containers.Map();
pdfPrefixMap(conditionPairs{1}{1}) = 'SAT_ERROR_CHOICE_';
pdfPrefixMap(conditionPairs{2}{1}) = 'SAT_ERROR_TIMING_';

cellPairInfo = load(spkCountFile,'cellPairInfo');
cellPairInfo = cellPairInfo.cellPairInfo;
spikeCorr = load(spkCountFile,'spikeCorr');
spikeCorr = spikeCorr.spikeCorr;

    %% compute the min-max for axis scaling
    % common for all psth plots
    allRast = [spikeCorr.xRasters;spikeCorr.yRasters];
    maxSpkPerSec = max(cellfun(@(x) max(fx_vecSmooth(mean(x),smoothBinWidthMs)),allRast));
    maxSpkPerSec = round(maxSpkPerSec*1000);
    psthYLims = [0 maxSpkPerSec];
    psthYTicks = [0 maxSpkPerSec];
    psthYTickLabel =  arrayfun(@(x) num2str(x,'%d'),psthYTicks','UniformOutput',false);
    psthYTickLabel(psthYTicks==0) = {'0'};
    psthYaxisLabel = 'Spk/s';
    % common for all Rsc plots
    rho_pval_idx = find(~cellfun(@isempty,regexp(spikeCorr.Properties.VariableNames,'rho_pval*','match')));
    allRsc = {};
    allRsc = arrayfun(@(x) [allRsc{:};spikeCorr{:,x}],rho_pval_idx(:),'UniformOutput',false);
    allRsc = vertcat(allRsc{:});
    minMaxRsc = cell2mat(cellfun(@(x) minmax(x(:,1)'),allRsc,'UniformOutput',false));
    minMaxRsc(minMaxRsc<-1 | minMaxRsc > 1) = NaN;
    minMaxRsc = minmax(minMaxRsc(:)');
    rscYlims(1) = round(minMaxRsc(1),1);
    rscYlims(2) = round(minMaxRsc(2),1);
    rscYTicks = [rscYlims(1) 0 rscYlims(2)];
    rscYTickLabel =  arrayfun(@(x) num2str(x,'%0.1f'),rscYTicks','UniformOutput',false);
    rscYTickLabel(rscYTicks==0) = {'0'};
    rscYaxisLabel = 'r_{sc}';
    
    
%%
[~,pdfBaseFile] = fileparts(spkCountFile);
colormap('jet');close gcf;
for cc = 1:numel(conditionPairs)
    
    conditions = conditionPairs{cc};  
    outPdfFile = fullfile(pdfOutputDir, [pdfPrefixMap(conditions{1}) pdfBaseFile '.pdf']);    
    %% plot each pair of conditions
    parentFig = getFigHandle();
    H_out = struct();
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);
    offsetsX=[0.004 (1:3).*0.248]; % for 4 columns
    offsetsY = [0.88 0.43]; %[0.90 0.45]; % for 2 rows
    startPos = 0.015; % top position of yPsth
    psthH = 0.05; psthW = psthH*3.5;
    gutter = 0.008; % space between plots
    
    %% plot each condition in a row
    for rowNum = 1:2
        condition = conditions{rowNum};
        axColor = [0.5 0.5 0.5];
        alignNames = unique(spikeCorr.alignedName,'stable');       
        for colNum = 1:3
            currSpkCorr = spikeCorr(strcmp(spikeCorr.condition,condition) ...
                & strcmp(spikeCorr.alignedName,alignNames{colNum}),:);
            unitSumm ={
                sprintf('%10s %10s %10s','Unit','nTrials','nSpikes')
                sprintf('%10s %10d %10d','X-Unit',size(currSpkCorr.xRasters{1},1),sum(currSpkCorr.xRasters{1}(:)))
                sprintf('%10s %10d %10d','Y-Unit',size(currSpkCorr.yRasters{1},1),sum(currSpkCorr.yRasters{1}(:)))
                };
            % bins for xUnit , yUnit
            psthBins = currSpkCorr.rasterBins{1};
            psthXLims = [min(psthBins) max(psthBins)];
            psthXTicks = min(psthBins):100:max(psthBins);
            psthXTickLabel =  arrayfun(@(x) num2str(x/1000,'%.1f'),psthXTicks','UniformOutput',false);
            psthXTickLabel(2:end-1) = {''};
            psthXTickLabel(psthXTicks==0) = {'0'};
            psthXaxisLabel = ['Time from ', currSpkCorr.alignedEvent{1},' (s)'];
            psthXaxisLabel = strrep(psthXaxisLabel,'SaccadePrimaryTempo','Saccade');
            
            %% H_Psth1
            pos(1) = offsetsX(colNum) + startPos;
            pos(2) = offsetsY(rowNum);
            pos(3:4) = [psthW psthH];
            H_out.H_psth1=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth1');
            
            rasters = currSpkCorr.xRasters{1};
            plot(psthBins,fx_vecSmooth(mean(rasters)*1000,smoothBinWidthMs),'LineWidth',1.5);
            annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
            doYLabel(gca,{'X-Unit'; psthYaxisLabel})
            hold on
            PlotUtils.plotRasters(rasters,psthBins);
            %% H_Psth2
            %pos(1) = offsetsX(colNum) + startPos;
            pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter);
            pos(3:4) = [psthW psthH];
            H_out.H_psth2=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth2');
            
            rasters = currSpkCorr.yRasters{1};
            plot(psthBins,fx_vecSmooth(mean(rasters)*1000,smoothBinWidthMs),'LineWidth',1.5);
            annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
            doYLabel(gca,{'Y-Unit'; psthYaxisLabel})
            hold on
            PlotUtils.plotRasters(rasters,psthBins);
            %% H_rsc50
            %pos(1) = offsetsX(colNum) + startPos;
            pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*2;
            pos(3:4) = [psthW psthH];
            H_out.H_psth2=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth2');
            rho_pval = currSpkCorr.rho_pval_50ms{1};
            sig_05 = currSpkCorr.critRho05;
            plot(psthBins,fx_vecSmooth(rho_pval(:,1),smoothBinWidthMs),'LineWidth',1.5);
            annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
            doYLabel(gca,'r_{sc} 50ms')
            hold on
            line(get(gca,'XLim'),[0 0],'Color','k')
            line(get(gca,'XLim'),[sig_05 sig_05],'Color','r')
            line(get(gca,'XLim'),[-sig_05 -sig_05],'Color','r')
            %sigIdx = rho_pval(:,2)<=0.05;
            %plot(psthBins(sigIdx),rho_pval(sigIdx,1),'r','LineWidth',2.0);
            
            %% H_rsc100
            %pos(1) = offsetsX(colNum) + startPos;
            pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*3;
            pos(3:4) = [psthW psthH];
            H_out.H_psth2=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth2');
            rho_pval = currSpkCorr.rho_pval_100ms{1};
            sig_05 = currSpkCorr.critRho05;
            plot(psthBins,fx_vecSmooth(rho_pval(:,1),smoothBinWidthMs),'LineWidth',1.5);
            annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
            doYLabel(gca,'r_{sc} 100ms')
            hold on
            line(get(gca,'XLim'),[0 0],'Color','k')
            line(get(gca,'XLim'),[sig_05 sig_05],'Color','r')
            line(get(gca,'XLim'),[-sig_05 -sig_05],'Color','r')
            %sigIdx = rho_pval(:,2)<=0.05;
            %plot(psthBins(sigIdx),rho_pval(sigIdx,1),'r','LineWidth',2.0);

            %% H_rsc200
            %pos(1) = offsetsX(colNum) + startPos;
            pos(2) = pos(2) - (psthH + gutter); % offsetsY(rowNum) - (psthH + gutter)*4;
            pos(3:4) = [psthW psthH];
            H_out.H_psth2=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth2');
            rho_pval = currSpkCorr.rho_pval_200ms{1};
            sig_05 = currSpkCorr.critRho05;
            plot(psthBins,fx_vecSmooth(rho_pval(:,1),smoothBinWidthMs),'LineWidth',1.5);
            annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
            doYLabel(gca,'r_{sc} 200ms')
            hold on
            line(get(gca,'XLim'),[0 0],'Color','k')
            line(get(gca,'XLim'),[sig_05 sig_05],'Color','r')
            line(get(gca,'XLim'),[-sig_05 -sig_05],'Color','r')
            %sigIdx = rho_pval(:,2)<=0.05;
            %plot(psthBins(sigIdx),rho_pval(sigIdx,1),'r','LineWidth',2.0);

            %% H_rsc499
            %pos(1) = offsetsX(colNum) + startPos;
            pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*5;
            pos(3:4) = [psthW psthH];
            H_out.H_psth2=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth2');
            rho_pval = currSpkCorr.rho_pval_400ms{1};
            sig_05 = currSpkCorr.critRho05;
            plot(psthBins,fx_vecSmooth(rho_pval(:,1),smoothBinWidthMs),'LineWidth',1.5);
            annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0,axColor);
            doYLabel(gca,'r_{sc} 400ms')
            hold on
            line(get(gca,'XLim'),[0 0],'Color','k')
            line(get(gca,'XLim'),[sig_05 sig_05],'Color','r')
            line(get(gca,'XLim'),[-sig_05 -sig_05],'Color','r')
            %sigIdx = rho_pval(:,2)<=0.05;
            %plot(psthBins(sigIdx),rho_pval(sigIdx,1),'r','LineWidth',2.0);
            doXLabel(gca,psthXaxisLabel);

            % Add unit summary annotation here
            annotation('textbox','Position',[pos(1) pos(2)-(psthH+gutter*2) 0.02 0.05],'String',char(unitSumm),...
                'FontSize',8,'FontWeight','bold','FitBoxToText','on','Interpreter','none','EdgeColor','none');

        end
    end
    addAnnotations(cellPairInfo.Pair_UID{1},outPdfFile, [cellPairInfo.X_area{1} ' vs ' cellPairInfo.Y_area{1}],...
                   conditions,alignNames);
    H_pairInfo = axes('parent',parentFig,'position',[0.01 0.01 0.98 0.06],...
        'box','on','XTick',[],'YTick',[],'layer','top','Tag','H_jpsthInfo');
    addPairInfo(H_pairInfo, cellPairInfo);
    drawnow;
    if savePdfFlag
        saveFigAs(outPdfFile);
        delete(parentFig);
    end
end

end

function addPairInfo(H_axes,cellPairInfo)
  cellInfo = cellPairInfo(1,contains(cellPairInfo.Properties.VariableNames,'X_'));
  cellInfo.Properties.VariableNames = strrep(cellInfo.Properties.VariableNames,'X_','');
  cellInfo(2,:) = cellPairInfo(1,contains(cellPairInfo.Properties.VariableNames,'Y_'));
  plotAddPairInfo(H_axes,cellInfo);
end

function addAnnotations(pairUid,pdfFile,xyAreas,rowNames,colNames)
% 
annotation('textbox',[0.10 0.95 0.05 0.05],'String',pairUid,'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
[~,fn,ext] = fileparts(pdfFile);
annotation('textbox',[0.35 0.95 0.05 0.05],'String',[fn ext],'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
annotation('textbox',[0.75 0.95 0.05 0.05],'String',xyAreas,'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
% conditions / alignNames
annotation('textbox',[0.01 0.93 0.05 0.05],'String',rowNames{1},'FontSize',20,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','green')
annotation('textbox',[0.12 0.91 0.05 0.05],'String',colNames{1},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
annotation('textbox',[0.32 0.91 0.05 0.05],'String',colNames{2},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
annotation('textbox',[0.55 0.91 0.05 0.05],'String',colNames{3},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
annotation('textbox',[0.01 0.48 0.05 0.05],'String',rowNames{2},'FontSize',20,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','red')

end
function doXLabel(H_axis,xLabel)
 yLim= get(H_axis,'Ylim');
 xPos = mean(get(H_axis,'Xlim')); % center in x
 if strcmp(get(H_axis,'YDir'),'reverse')
    yPos = yLim(2) + range(yLim)*0.05;
 else
    yPos = yLim(1) - range(yLim)*0.05; % 5% below the x-axis
 end
 
 xlabel(xLabel,'Position',[xPos yPos],'VerticalAlignment','top',...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold',...
     'FontAngle', 'italic','Color','black'); 
end

function doYLabel(H_axis,yLabel)
 xLim= get(H_axis,'Xlim');
 yPos = mean(get(H_axis,'Ylim')); % center in y
 if strcmp(get(H_axis,'YAxisLocation'),'right')
     xPos = xLim(2) + range(xLim)*0.02; % 10% to right
 else
     v = get(H_axis,'View');
     if (v(1) == 0) % normal view
         xPos = xLim(1) - range(xLim)*0.08; % left
     else
         xPos = xLim(1) - range(xLim)*0.01; % down
     end
 end
 ylabel(yLabel,'Position',[xPos yPos],'VerticalAlignment','top',...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold',...
     'FontAngle', 'italic','Color','black'); 
end


function annotateAxis(H_axis,xOry,lims_,ticks_,labels_,rotDeg,axColor)
if xOry=='y'
    set(H_axis,'ylim',lims_,'ytick',ticks_,'YTickLabel',labels_,...
        'YTickLabelRotation',rotDeg,'ycolor',[0.5 0.5 0.5],'YColor',axColor);
elseif xOry=='x'
    set(H_axis,'xlim',lims_,'xtick',ticks_,'XTickLabel',labels_,...
        'XTickLabelRotation',rotDeg,'xcolor',[0.5 0.5 0.5],'XColor',axColor);
    set(H_axis,'XGrid','on','GridLineStyle','--', 'GridColor', [0.3 0.3 0.3])
    line([0 0],lims_,'Color',[0.1 0.1 0.1],'LineStyle',':')
end
end

%%
function [H_Figure] = getFigHandle()
set(0, 'DefaultFigureColormap', jet(64));
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






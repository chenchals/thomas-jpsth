function [] = corrSpkCountPlot(spkCountFile,pdfOutputDir,savePdfFlag)
%CORRSPKCOUNTPLOT Summary of this function goes here

    %% Put this in the function....
    smoothBinWidthMs = 5;
    fx_vecSmooth = @(x,w) smoothdata(x,'movmean',w,'omitnan');
    fx_chanNo = @(x) str2double(char(regexp(x,'(^\d{1,2})','match')));
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

    %%
    [~,pdfBaseFile] = fileparts(spkCountFile);
    [alignNames,idx] = unique(spikeCorr.alignedName,'stable');
    rhoPvalWins = spikeCorr.rho_pval_win(idx);
    %colormap('jet');close gcf;
    colrs = lines;
    for cc = 1:numel(conditionPairs)

        conditions = conditionPairs{cc};
        outPdfFile = fullfile(pdfOutputDir, [pdfPrefixMap(conditions{1}) pdfBaseFile '.pdf']);
        %% plot each pair of conditions
        parentFig = getFigHandle();
        H_out = struct();
        ss = get(0,'ScreenSize');
        aspectRatio = ss(3)/ss(4);
        offsetsX=[0.01 (1:4).*0.17]; % for 5 columns
        offsetsY = [0.88 0.45]; %[0.90 0.45]; % for 2 rows
        startPos = 0.017; % top position of yPsth
        psthH = 0.05; psthW = psthH*2.5;
        gutter = 0.01; % space between plots

        %% plot each condition in a row
        for rowNum = 1:2
            condition = conditions{rowNum};
            axColor = [0.5 0.5 0.5];
            for colNum = 1:4
                rhoPvalWin = rhoPvalWins{colNum};
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
                plot(psthBins,fx_vecSmooth(mean(rasters)*1000,smoothBinWidthMs),'LineWidth',1);
                annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,0,axColor);
                annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
                doYLabel(gca,{'X-Unit'; psthYaxisLabel})
                hold on
                sortMarkers = getAlignedSortMarkers(currSpkCorr);
                PlotUtils.plotRasters(rasters,psthBins,sortMarkers);
                addPatch(gca,rhoPvalWin);
                
                %% H_Psth2
                %pos(1) = offsetsX(colNum) + startPos;
                pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter);
                pos(3:4) = [psthW psthH];
                H_out.H_psth2=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_psth2');

                rasters = currSpkCorr.yRasters{1};
                plot(psthBins,fx_vecSmooth(mean(rasters)*1000,smoothBinWidthMs),'LineWidth',1);
                annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,0,axColor);
                annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
                doYLabel(gca,{'Y-Unit'; psthYaxisLabel})
                hold on
                PlotUtils.plotRasters(rasters,psthBins,sortMarkers);
                addPatch(gca,rhoPvalWin);

                %% H_rsc50
                %pos(1) = offsetsX(colNum) + startPos;
                pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*2;
                pos(3:4) = [psthW psthH];
                H_out.H_rsc50=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_rsc50');

                rho_pval = currSpkCorr.rho_pval_50ms{1};
                rho_pvalZ = currSpkCorr.rho_pval_50ms_Z{1};
                plotRhoPvals(psthBins,rho_pval,rho_pvalZ,smoothBinWidthMs,fx_vecSmooth,colrs);

                annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
                annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
                doYLabel(gca,'r_{sc} 50ms')
                line(get(gca,'XLim'),[0 0],'Color','k')
                addPatch(gca,rhoPvalWin);

                %% H_rsc100
                %pos(1) = offsetsX(colNum) + startPos;
                pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*3;
                pos(3:4) = [psthW psthH];
                H_out.H_rsc100=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_rsc100');

                rho_pval = currSpkCorr.rho_pval_100ms{1};
                rho_pvalZ = currSpkCorr.rho_pval_100ms_Z{1};
                plotRhoPvals(psthBins,rho_pval,rho_pvalZ,smoothBinWidthMs,fx_vecSmooth,colrs);

                annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
                annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
                doYLabel(gca,'r_{sc} 100ms')
                hold on
                line(get(gca,'XLim'),[0 0],'Color','k')
                addPatch(gca,rhoPvalWin);

                %% H_rsc200
                %pos(1) = offsetsX(colNum) + startPos;
                pos(2) = pos(2) - (psthH + gutter); % offsetsY(rowNum) - (psthH + gutter)*4;
                pos(3:4) = [psthW psthH];
                H_out.H_rsc200=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_rsc200');

                rho_pval = currSpkCorr.rho_pval_200ms{1};
                rho_pvalZ = currSpkCorr.rho_pval_200ms_Z{1};
                plotRhoPvals(psthBins,rho_pval,rho_pvalZ,smoothBinWidthMs,fx_vecSmooth,colrs);

                annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
                annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
                doYLabel(gca,'r_{sc} 200ms')
                hold on
                line(get(gca,'XLim'),[0 0],'Color','k')
                addPatch(gca,rhoPvalWin);

                %% H_rsc400
                %pos(1) = offsetsX(colNum) + startPos;
                pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*5;
                pos(3:4) = [psthW psthH];
                H_out.H_rsc400=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_rsc400');

                rho_pval = currSpkCorr.rho_pval_400ms{1};
                rho_pvalZ = currSpkCorr.rho_pval_400ms_Z{1};
                plotRhoPvals(psthBins,rho_pval,rho_pvalZ,smoothBinWidthMs,fx_vecSmooth,colrs);

                annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
                annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0,axColor);
                doYLabel(gca,'r_{sc} 400ms')
                hold on
                line(get(gca,'XLim'),[0 0],'Color','k')
                addPatch(gca,rhoPvalWin);

                doXLabel(gca,psthXaxisLabel);
                % Add unit summary annotation here
                annotation('textbox','Position',[pos(1) pos(2)-(psthH+gutter*2) 0.02 0.05],'String',char(unitSumm),...
                    'FontSize',8,'FontWeight','bold','FitBoxToText','on','Interpreter','none','EdgeColor','none');

            end
            %% Draw static window spike count corr
            staticCols  = {'xSpkCount_win','ySpkCount_win','rho_pval_win','rho_pval_static'};
            statColrIdx = [6 2];
            for z = 1:2
                if z > 1
                    staticCols = strcat(staticCols,'_Z');
                end
                scatColr = colrs(statColrIdx(z),:);
                % get min-max of spk counts for scaling
                maxSpkCountX = max(cell2mat(spikeCorr.(staticCols{1})));
                maxSpkCountY = max(cell2mat(spikeCorr.(staticCols{2})));
                if z == 1
                    maxSpkCountX = maxSpkCountX + mod(maxSpkCountX,2);
                    maxSpkCountY = maxSpkCountY + mod(maxSpkCountY,2);
                else
                    maxSpkCountX = round(maxSpkCountX,2);
                    maxSpkCountY = round(maxSpkCountY,2);
                end
                spkCountLimsX = [0 maxSpkCountX];
                spkCountLimsY = [0 maxSpkCountY];
                spkCountTicksX = 0:maxSpkCountX/4:maxSpkCountX;
                spkCountTicksY = 0:maxSpkCountY/4:maxSpkCountY;
                if z == 1
                    spkCountTickLabelX =  arrayfun(@(x) num2str(x,'%d'),spkCountTicksX','UniformOutput',false);
                    spkCountTickLabelY =  arrayfun(@(x) num2str(x,'%d'),spkCountTicksY','UniformOutput',false);
                    spkCountXAxisLabel = {'X-Unit', 'Spk.Count'};
                    spkCountYAxisLabel = {'Y-Unit', 'Spk.Count'};
                else
                    spkCountTickLabelX =  arrayfun(@(x) num2str(x,'%0.1f'),spkCountTicksX','UniformOutput',false);
                    spkCountTickLabelY =  arrayfun(@(x) num2str(x,'%0.1f'),spkCountTicksY','UniformOutput',false);
                    spkCountXAxisLabel = {'X-Unit', 'sum(z-score)'};
                    spkCountYAxisLabel = {'Y-Unit', 'sum(z-score)'};
                end
                spkCountTickLabelX(2:end-1) = repmat({' '},numel(spkCountTickLabelX)-2,1);
                spkCountTickLabelY(2:end-1) = repmat({' '},numel(spkCountTickLabelY)-2,1);

                pltW = psthW/2.6;
                pltH = pltW*aspectRatio;

                for s = 1:4
                    % Draw 3 scatter plots across for Baseline, visual, postsac
                    pos(1) = offsetsX(5) + gutter + (pltW + gutter*3)*(s-1);
                    pos(2) = offsetsY(rowNum) - (psthH + gutter) - (pltH + gutter*6)*(z-1);
                    pos(3:4) = [pltW pltH];
                    H_out.H_rscBl=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_rscBl');
                    plotData = spikeCorr(strcmp(spikeCorr.condition,condition) ...,
                        & strcmp(spikeCorr.alignedName,alignNames{s}),...
                        staticCols);
                    x = plotData.(staticCols{1}){1};
                    y = plotData.(staticCols{2}){1};
                    scatter(x,y,20,repmat(scatColr,numel(x),1),'filled','o');
                    hold on
                    annotateAxis(gca,'y',spkCountLimsY,spkCountTicksY,spkCountTickLabelY,0,axColor);
                    annotateAxis(gca,'x',spkCountLimsX,spkCountTicksX,spkCountTickLabelX,0,axColor);
                    set(gca,'YGrid','on','GridLineStyle','--', 'GridColor', [0.3 0.3 0.3])

                    doYLabel(gca,spkCountYAxisLabel)
                    doXLabel(gca,spkCountXAxisLabel)
                    if z == 1
                        titleStr = {alignNames{s}...
                            ['[' num2str(plotData.(staticCols{3}){1},'%d ') ']' ]...
                            sprintf('\\rho = %0.2f, p = %0.2e', ...
                            plotData.(staticCols{4}){1}(1),plotData.(staticCols{4}){1}(2))};
                    else
                        titleStr = {...
                            sprintf('\\rho = %0.2f, p = %0.2e', ...
                            plotData.(staticCols{4}){1}(1),plotData.(staticCols{4}){1}(2))};
                    end
                    title(titleStr,'Interpreter','tex')
                end
            end
            
            %% for each condition get waveforms and wave Widths for unitx and unit y
            xyWf = spikeCorr(strcmp(spikeCorr.condition,condition),{'xWaves','yWaves'});
            unitsTxt = {cellPairInfo.X_unit{1} cellPairInfo.Y_unit{1}};
            wColrs = {'r','c'};
            % at 40Kz = 25 microsec/sample
            binSize = (1/40);
            wfBins = (0:size(xyWf.xWaves{1}{1},2))*binSize; % in millisecs
            wfXlim = [min(wfBins) max(wfBins)];
            wfXTicks = min(wfBins):0.1:max(wfBins);
            wfXTickLabel =  arrayfun(@(x) num2str(x,'%0.2f'),wfXTicks','UniformOutput',false);
            wfXTickLabel(wfXTicks==0) = {'0'};
            wfXaxisLabel = 'Time (ms)';          
            wfYaxisLabel = 'AD Units';
             
            %% for each condition plot aggregated waveform for unitx and unit y
            pos(1) = offsetsX(5) + gutter;
            pos(2) = offsetsY(rowNum) - (pltH + gutter*6)*2.5;
            pos(3:4) = [pltW*2.9 pltH*1.3];
            H_out.H_wav=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_wav');
            plotWaveforms(cell2mat(xyWf.xWaves{1}),binSize,0,wColrs{1});
            hold on
            plotWaveforms(cell2mat(xyWf.yWaves{1}),binSize,0,wColrs{2});
            doYLabel(gca,wfYaxisLabel)          
            annotateAxis(gca,'x',wfXlim,wfXTicks,wfXTickLabel,0,axColor);
            doXLabel(gca,wfXaxisLabel)            
            legend(unitsTxt,'Location', 'northwest','Box','off')
           
        %% end rowNum
        end
        chanStr = 'Different Channels';
        if fx_chanNo(cellPairInfo.X_unit{1}) == fx_chanNo(cellPairInfo.Y_unit{1})
            chanStr = 'Same Channel';
        end
        addAnnotations(cellPairInfo.Pair_UID{1},outPdfFile, [cellPairInfo.X_area{1} ' vs ' cellPairInfo.Y_area{1}],...
            chanStr,conditions,alignNames);
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


function [] = plotRhoPvals(psthBins,rho_pval,rho_pvalZ,smoothBinWidthMs,fx_handle,colrs)
    yVals = fx_handle(rho_pval(:,1),smoothBinWidthMs);
    p = plot(psthBins,yVals,'Color',colrs(6,:));
    hold on
    yVals(rho_pval(:,2)>0.05) = NaN;
    plot(psthBins,yVals,'Color',p.Color,'LineWidth',5);
    yValsZ = fx_handle(rho_pvalZ(:,1),smoothBinWidthMs);
    p = plot(psthBins,yValsZ,'Color',colrs(2,:));
    yValsZ(rho_pvalZ(:,2)>0.05) = NaN;
    plot(psthBins,yValsZ,'Color',p.Color ,'LineWidth',5);
end

function [] = addPatch(H_ax,winMs)
  x = [winMs ;winMs];
  x=x(:)';
  y = get(H_ax,'ylim');
  y = [y fliplr(y)];
  c = [0.8 0.8 0.8];
  fill(x, y,c,'FaceAlpha',0.4,'LineStyle','none')
end

function addPairInfo(H_axes,cellPairInfo)
  cellInfo = cellPairInfo(1,contains(cellPairInfo.Properties.VariableNames,'X_'));
  cellInfo.Properties.VariableNames = strrep(cellInfo.Properties.VariableNames,'X_','');
  cellInfo(2,:) = cellPairInfo(1,contains(cellPairInfo.Properties.VariableNames,'Y_'));
  plotAddPairInfo(H_axes,cellInfo);
end

function addAnnotations(pairUid,pdfFile,xyAreas,chanStr,rowNames,colNames)
    %
    annotation('textbox',[0.10 0.95 0.05 0.05],'String',pairUid,'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    [~,fn,ext] = fileparts(pdfFile);
    annotation('textbox',[0.35 0.95 0.05 0.05],'String',[fn ext],'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    annotation('textbox',[0.70 0.95 0.05 0.05],'String',xyAreas,'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    h = annotation('textbox',[0.85 0.95 0.10 0.05],'String',chanStr,'FontSize',24,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none');
    if contains(chanStr,'Same')
        set(h,'Color','r');
    end% conditions / alignNames
    annotation('textbox',[0.01 0.93 0.05 0.05],'String',rowNames{1},'FontSize',20,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','green')
    annotation('textbox',[0.07 0.91 0.05 0.05],'String',colNames{1},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.22 0.91 0.05 0.05],'String',colNames{2},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.38 0.91 0.05 0.05],'String',colNames{3},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.56 0.91 0.05 0.05],'String',colNames{4},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
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
         xPos = xLim(1) - range(xLim)*0.05; % left
     else
         xPos = xLim(1) - range(xLim)*0.01; % down
     end
 end
 ylabel(yLabel,'Position',[xPos yPos],'VerticalAlignment','bottom',...
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

function [sortMarkers] = getAlignedSortMarkers(thisTbl)
   % CueOn time is always 3500.
    fx_getFirstSortColIdx = @(tbl) find(strcmp(tbl.Properties.VariableNames,'firstSortBy'));
    fx_getSecondSortColIdx = @(tbl) find(strcmp(tbl.Properties.VariableNames,'secondSortBy'));
    alignTimeMinusCueOnTime = thisTbl.alignTime{1} - 3500;
    
    sortMarkers = cell(2,1);
    if ~isempty(fx_getFirstSortColIdx(thisTbl))
        sortMarkers{1} = thisTbl{1,fx_getFirstSortColIdx(thisTbl)}{1} - alignTimeMinusCueOnTime;
    end
    if ~isempty(fx_getSecondSortColIdx(thisTbl))
        sortMarkers{2} = thisTbl{1,fx_getSecondSortColIdx(thisTbl)}{1} - alignTimeMinusCueOnTime;
    end
end





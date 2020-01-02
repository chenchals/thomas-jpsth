function [] = corrSpkCountPlot(spkCountFile,pdfOutputDir,savePdfFlag)
%CORRSPKCOUNTPLOT Summary of this function goes here
    %% ....
    fs = 6;
    if ismac
        fs = 8;
    end
    plot400MsMovingWin = false;
    doErrorRewardGrade = true;
    smoothBinWidthMs = 5;
    fx_vecSmooth = @(x,w) smoothdata(x,'movmean',w,'omitnan');
    conditionPairs = {
        {'FastErrorChoice','AccurateErrorChoice'};
        {'FastErrorTiming','AccurateErrorTiming'};
        {'FastCorrect','AccurateCorrect'};
        };
    pdfPrefixMap = containers.Map();
    pdfPrefixMap(conditionPairs{1}{1}) = 'SAT_ERROR_CHOICE_';
    pdfPrefixMap(conditionPairs{2}{1}) = 'SAT_ERROR_TIMING_';
    pdfPrefixMap(conditionPairs{3}{1}) = 'SAT_CORRECT_';

    cellPairInfo = load(spkCountFile,'cellPairInfo');
    cellPairInfo = cellPairInfo.cellPairInfo;
    spikeCorr = load(spkCountFile,'spkCorr');
    spikeCorr = spikeCorr.spkCorr;

    %% compute the min-max for axis scaling
    % common for all psth plots
    % choose static time win to use: 50|150|200
    staticMsToUse = 150;
    allRast = [spikeCorr.xRasters;spikeCorr.yRasters];
    maxSpkPerSec = max(cellfun(@(x) max(fx_vecSmooth(mean(x),smoothBinWidthMs)),allRast));
    maxSpkPerSec = round(maxSpkPerSec*1000);
    psthYLims = [0 maxSpkPerSec];
    psthYTicks = [0 maxSpkPerSec];
    psthYTickLabel =  arrayfun(@(x) num2str(x,'%d'),psthYTicks','UniformOutput',false);
    psthYTickLabel(psthYTicks==0) = {'0'};
    psthYaxisLabel = 'Spk. Count';
    % common for all Rsc plots
    rho_pval_idx = find(~cellfun(@isempty,regexp(spikeCorr.Properties.VariableNames,'rho_pval_(50|100|200)ms','match')));
    rho_pval_cols = spikeCorr.Properties.VariableNames(rho_pval_idx)';
    allRsc = {};
    allRsc = arrayfun(@(x) [allRsc{:};spikeCorr{:,x}],rho_pval_idx(:),'UniformOutput',false);
    allRsc = vertcat(allRsc{:});
    minMaxRsc = cell2mat(cellfun(@(x) minmax(x(:,1)'),allRsc,'UniformOutput',false));
    minMaxRsc(minMaxRsc<-1 | minMaxRsc > 1) = NaN;
    minMaxRsc = minmax(minMaxRsc(:)');
    rscYlims(1) = round(minMaxRsc(1)-0.05,1);
    rscYlims(2) = round(minMaxRsc(2)+0.05,1);
    rscYTicks = [rscYlims(1) 0 rscYlims(2)];
    rscYTickLabel =  arrayfun(@(x) num2str(x,'%0.1f'),rscYTicks','UniformOutput',false);
    rscYTickLabel(rscYTicks==0) = {'0'};

    %%
    rhoPvalStaticWinCol = strcat('rho_pval_win',num2str(staticMsToUse,'_%dms'));
    [~,pdfBaseFile] = fileparts(spkCountFile);
    [alignNames,idx] = unique(spikeCorr.alignedName,'stable');
    rhoPvalWins = spikeCorr.(rhoPvalStaticWinCol)(idx);
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
        offsetsY = [0.86 0.42]; %[0.90 0.45]; % for 2 rows
        startPos = 0.017; % top position of yPsth
        % For 4 timewins 
        if plot400MsMovingWin
          psthH = 0.05; psthW = psthH*2.5; %#ok<*UNRCH>
        else
          psthH = 0.06; psthW = 0.125;
        end
        gutter = 0.01; % space between plots

        %% plot each condition in a row
        for rowNum = 1:2
            condition = conditions{rowNum};
            axColor = [0.5 0.5 0.5];
            % it is possible that a particular condition may not exist
            % example FastErrorTiming for SEF-SC --> PAIR_0590
            if sum(strcmp(spikeCorr.condition,condition)) == 0
                continue;
            end
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
                %psthXaxisLabel = strrep(psthXaxisLabel,'SaccadePrimaryTempo','Saccade');

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
                
                %% Plot rho-spike counts Z-scored by trial
                for n = 1:3 % 3 plots [50ms, 100ms,200ms]
                    pos(2) = pos(2) - (psthH + gutter); %offsetsY(rowNum) - (psthH + gutter)*2;
                    pos(3:4) = [psthW psthH];
                    colName = rho_pval_cols{n};
                    suff = regexp(colName,'\d*ms','match');
                    tag = strcat('H_rsc_',suff{1});
                    ylabelTxt = {['\fontsize{12} \rho \fontsize{' num2str(fs) '} Spk.Corr.'], suff{1}};
                    H_out.(tag) = axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag',tag);
                    rho_pval = currSpkCorr.(colName){1};
                    plotRhoPvals(psthBins,rho_pval,[],smoothBinWidthMs,fx_vecSmooth,colrs);
                    annotateAxis(gca,'y',rscYlims,rscYTicks,rscYTickLabel,0,axColor);
                    annotateAxis(gca,'x',psthXLims,psthXTicks,{},0,axColor);
                    doYLabel(gca,ylabelTxt)
                    line(get(gca,'XLim'),[0 0],'Color','k')
                    addPatch(gca,rhoPvalWin);
                    
                end
                doXLabel(gca,psthXaxisLabel);
                
                set(gca,'XTickLabel',psthXTickLabel);
                
                % Add unit summary annotation here
                annotation('textbox','Position',[pos(1) pos(2)-(psthH*0.9) 0.02 0.04],'String',char(unitSumm),...
                    'FontSize',fs,'FontWeight','bold','FitBoxToText','on','Interpreter','none','EdgeColor','none');
                
            end
            %% Draw static window spike count corr
            staticCols  = {'xSpkCount_win','ySpkCount_win','rho_pval_win','rho_pval_static'};
            staticCols = strcat(staticCols,num2str(staticMsToUse,'_%dms'));
            currSpkCorr = spikeCorr(strcmp(spikeCorr.condition,condition),:);
            statColrIdx = [6 2];
            for z = 1:1
                if z > 1
                    staticCols = strcat(staticCols,'_Z_trial');
                end
                scatColr = colrs(statColrIdx(z),:);
                % get min-max of spk counts for scaling
                maxSpkCountX = max(cell2mat(currSpkCorr.(staticCols{1})));
                maxSpkCountY = max(cell2mat(currSpkCorr.(staticCols{2})));
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
                    spkCountXAxisLabel = {'X-Spk.Count'};
                    spkCountYAxisLabel = {'Y-Spk.Count'};
                else
                    spkCountTickLabelX =  arrayfun(@(x) num2str(x,'%0.1f'),spkCountTicksX','UniformOutput',false);
                    spkCountTickLabelY =  arrayfun(@(x) num2str(x,'%0.1f'),spkCountTicksY','UniformOutput',false);
                    spkCountXAxisLabel = {'X-sum(z-score)'};
                    spkCountYAxisLabel = {'Y-sum(z-score)'};
                end
                spkCountTickLabelX(2:end-1) = repmat({' '},numel(spkCountTickLabelX)-2,1);
                spkCountTickLabelY(2:end-1) = repmat({' '},numel(spkCountTickLabelY)-2,1);

                pltW = psthW/2.6;
                pltH = pltW*aspectRatio;

                for s = 1:4
                    % Draw scatter plots across for Baseline, visual,
                    % postsac, postRew
                    pos(1) = offsetsX(5) + gutter + (pltW + gutter*3)*(s-1);
                    pos(2) = offsetsY(rowNum) - (psthH*.5) - (pltH + gutter*4)*(z-1) - pltH*0.5;
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
                    
                    if plotData.(staticCols{4}){1}(2) <= 0.05
                        signif = '^*';
                        titleColor = 'r';
                    else
                        signif = '';
                        titleColor = 'k';
                    end                  
                    if z == 1
                        titleStr = {alignNames{s}...
                            ['[' num2str(plotData.(staticCols{3}){1},'%d ') ']' ]...
                             sprintf(['\\rho = %0.2f, p = %0.2e' signif], ...
                            plotData.(staticCols{4}){1}(1),plotData.(staticCols{4}){1}(2))};
                    else
                        titleStr = {...
                            sprintf(['\\rho = %0.2f, p = %0.2e' signif], ...
                            plotData.(staticCols{4}){1}(1),plotData.(staticCols{4}){1}(2))};
                    end
                    title(titleStr,'Interpreter','tex','Color',titleColor)
                end
            end
            
            %% for each condition get waveforms for unitx and unit y
            xyWf = spikeCorr(strcmp(spikeCorr.condition,condition),{'xWaves','yWaves'});
            %wColrs = {'r','c'};
            wColrs = {[1 0 0],[0 1 1]};
            % at 40Kz = 25 microsec/sample
            binSize = (1/40);
            wfBins = (0:size(xyWf.xWaves{1}{1},2))*binSize; % in millisecs
            wfXlim = [min(wfBins) max(wfBins)];
            wfXTicks = min(wfBins):0.1:max(wfBins);
            wfXTickLabel =  arrayfun(@(x) num2str(x,'%0.2f'),wfXTicks','UniformOutput',false);
            %wfXTickLabel(wfXTicks==0) = {'0'};
            wfXTickLabel(wfXTicks==0) = {'0'};
            wfXTickLabel(2:end-1) = {' '};
            
            wfXaxisLabel = 'Time (ms)';          
            wfYaxisLabel = 'AD Units';
            % waveforms across all alignments
            xWaves = cell2mat(cellfun(@(x) cell2mat(x),xyWf.xWaves,'UniformOutput',false));
            yWaves = cell2mat(cellfun(@(x) cell2mat(x),xyWf.yWaves,'UniformOutput',false));
            unitsTxt = {sprintf('%s (%d)',cellPairInfo.X_unit{1},size(xWaves,1)),...
                sprintf('%s (%d)',cellPairInfo.Y_unit{1},size(yWaves,1))};
            % for each condition plot aggregated waveform for unitx and unit y
            pos(1) = offsetsX(5) + gutter;
            pos(2) = offsetsY(rowNum) - (pltH + gutter*4)*2.5 + pltH*0.5;

            pos(3:4) = [pltW*1.6 pltH*1.3];
            H_out.H_wav=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_wav');
            plotWaveforms(xWaves,binSize,0,wColrs{1});
            hold on
            plotWaveforms(yWaves,binSize,0,wColrs{2});
            doYLabel(gca,wfYaxisLabel)          
            annotateAxis(gca,'x',wfXlim,wfXTicks,wfXTickLabel,0,axColor);
            doXLabel(gca,wfXaxisLabel)            
            legend(unitsTxt,'Location', 'northwest','Box','off')
            
            %% for each condition get waveforms for unitx and unit y
            sampleTime = 1000/40; % 40KHz, = 25 microSecs
            xyWavWidths = spikeCorr(strcmp(spikeCorr.condition,condition),{'xWaveWidths','yWaveWidths'});
            xWavWidths = cell2mat(cellfun(@(x) cell2mat(x),xyWavWidths.xWaveWidths,'UniformOutput',false));
            yWavWidths = cell2mat(cellfun(@(x) cell2mat(x),xyWavWidths.yWaveWidths,'UniformOutput',false));
            xWavWidths = abs(xWavWidths)*sampleTime;
            yWavWidths = abs(yWavWidths)*sampleTime;
            wWXlim = [0 max([max(xWavWidths),max(yWavWidths)])+100];
            wWXTicks = min(wWXlim):50:max(wWXlim);
            wWXTickLabel =  arrayfun(@(x) num2str(x,'%d'),wWXTicks','UniformOutput',false);
            wWXTickLabel(contains(wWXTickLabel,'50')) = {' '};
            wWXTickLabel(wWXTicks==0) = {'0'};
            wWXTickLabel(3:end-2) = {' '};
            wWXaxisLabel = 'Time (microsec)';          
            wWYaxisLabel = 'Frequency';
            unitsTxt = {cellPairInfo.X_unit{1},cellPairInfo.Y_unit{1}};
  
            pos(1) = pos(1) + pos(3) + gutter*3;
            pos(2) = pos(2);
            pos(3:4) = [pltW*1.6 pltH*1.3];
            H_out.H_wavWid=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_wavWid');
            
            histogram(xWavWidths,50,'BinLimits',wWXlim,'FaceColor',wColrs{1});
            hold on
            histogram(yWavWidths,50,'BinLimits',wWXlim,'FaceColor',wColrs{2});
            hold off
            doYLabel(gca,wWYaxisLabel)
            annotateAxis(gca,'x',wWXlim,wWXTicks,wWXTickLabel,0,axColor);
            doXLabel(gca,wWXaxisLabel)
            legend(unitsTxt,'Location', 'northwest','Box','off')
            % manage tick labels
            ticklabels = get(gca,'YTickLabel');
            ticklabels(2:end-1) = {' '};
            set(gca,'YTickLabel',ticklabels);
            
            % Mean std of spike widths
            xMean = mean(xWavWidths);
            xStd = std(xWavWidths);
            yMean = mean(yWavWidths);
            yStd = std(yWavWidths);
            unitsTxt = categorical({cellPairInfo.X_unit{1},cellPairInfo.Y_unit{1}});
            yLims = [0 max([xMean+xStd, yMean+yStd]+100)]; 

            pos(1) = pos(1) + pos(3) + gutter*3;
            pos(2) = pos(2);
            pos(3:4) = [pltW*1.6 pltH*1.3];
            H_out.H_wavWidM=axes('parent',parentFig,'position',pos,'box','on', 'layer','top','Tag','H_wavWidM');
            bar(unitsTxt(1),xMean,'FaceColor',wColrs{1});
            hold on
            h_ex = errorbar(unitsTxt(1),xMean,xStd,'.','Color',[0 0 0]);
            set(h_ex,'HandleVisibility','off');
            bar(unitsTxt(2),yMean,'FaceColor',wColrs{2});
            h_ey = errorbar(unitsTxt(2),yMean,yStd,'.','Color',[0 0 0]);
            set(h_ey,'HandleVisibility','off');
            
            ylim(yLims);
            ylabel('SpkeWidth (microsec.)','VerticalAlignment','bottom',...
                'HorizontalAlignment','center','FontSize',fs,'FontWeight','bold',...
                'FontAngle', 'italic','Color','black');
            
            meanStd = {sprintf('\\mu%4.1f\\pm %4.1f',xMean,xStd), sprintf('\\mu%4.f\\pm%4.1f',yMean,yStd)};
            %legend(meanStd,'Location','northeast','Box','off','FontSize',5)
            xTicks = get(gca,'XTick');
            yLocs = [xMean+1.2*xStd, yMean+1.2*yStd];
            text(xTicks,yLocs,meanStd,'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',fs,'FontWeight','bold'); 
            
            ticklabels = get(gca,'YTickLabel');
            ticklabels(2:end-1) = {' '};
            set(gca,'YTickLabel',ticklabels);
                        
        %% end rowNum
        end
        if cellPairInfo.isOnSameChannel
            chanStr = 'Same (dist: 0 mm)';
        else
            chanStr = sprintf('Different (dist: %0.3f mm)',cellPairInfo.XY_Dist{1});
        end
        addAnnotations(cellPairInfo.Pair_UID{1},outPdfFile, [cellPairInfo.X_area{1} ' vs ' cellPairInfo.Y_area{1}],...
            chanStr,conditions,alignNames);
        if doErrorRewardGrade
            annotateErrorRewardGrade(cellPairInfo);
        end
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

%%
function [] = annotateErrorRewardGrade(cellPairInfo)
        fs = 8;
        if ismac
            fs = 10;
        end
        errorRewardlabel1 ={
            sprintf('%15s ','Unit')
            sprintf('%15s ','\color{black} X-Unit')
            sprintf('%15s ','\color{black} Y-Unit')
            };
        errorRewardlabel2 ={
            sprintf('%15s ','IsErrorGrade')
            sprintf('%15s ',isTFTxt(cellPairInfo{1,'X_isErrGrade'}))
            sprintf('%15s ',isTFTxt(cellPairInfo{1,'Y_isErrGrade'}))
            };
        errorRewardlabel3 ={
            sprintf('%15s ','IsRewardGrade')
            sprintf('%15s ',isTFTxt(cellPairInfo{1,'X_isRewGrade'}))
            sprintf('%15s ',isTFTxt(cellPairInfo{1,'Y_isRewGrade'}))
            };
        % Add errorRewardlabel summary annotation here
        annotation('textbox','Position',[0.80 0.96 0.10 0.04],'String',char(errorRewardlabel1),...
            'Interpreter','tex','EdgeColor','none','FontWeight','bold','FitBoxToText','on','FontSize', fs);
        annotation('textbox','Position',[0.85 0.96 0.10 0.04],'String',char(errorRewardlabel2),...
            'Interpreter','tex','EdgeColor','none','FontWeight','bold','FitBoxToText','on', 'FontSize', fs);
        annotation('textbox','Position',[0.90 0.96 0.10 0.04],'String',char(errorRewardlabel3),...
            'Interpreter','tex','EdgeColor','none','FontWeight','bold','FitBoxToText','on','FontSize', fs);

end

function [txt] = isTFTxt(flag)
    if flag
        if ispc || ismac
            txt ='\color{green} YES';
        elseif ismac
            txt = ['\fontname{wingdings} \fontsize{30} \color{green}' char(252)];
        end
    else
        if ispc || ismac
            txt ='\color{red} NO';
        elseif ismac
            txt = ['\fontname{wingdings} \fontsize{30} \color{red}' char(251)];
        end
    end
end

function [] = plotRhoPvals(psthBins,rho_pval,rho_pvalZ,smoothBinWidthMs,fx_handle,colrs)
    if ~isempty(rho_pval)
        yVals = fx_handle(rho_pval(:,1),smoothBinWidthMs);
        p = plot(psthBins,yVals,'Color',colrs(6,:));
        hold on
        yVals(rho_pval(:,2)>0.05) = NaN;
        plot(psthBins,yVals,'Color',p.Color,'LineWidth',5);
    end
    if ~isempty(rho_pvalZ)
        yValsZ = fx_handle(rho_pvalZ(:,1),smoothBinWidthMs);
        hold on
        p = plot(psthBins,yValsZ,'Color',colrs(2,:));
        yValsZ(rho_pvalZ(:,2)>0.05) = NaN;
        plot(psthBins,yValsZ,'Color',p.Color ,'LineWidth',5);
    end
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
  cellInfo.XY_Dist = cellPairInfo.XY_Dist;
  cellInfo.Properties.VariableNames = strrep(cellInfo.Properties.VariableNames,'X_','');
  cellInfo(2,:) = cellPairInfo(1,contains(cellPairInfo.Properties.VariableNames,'Y_'));
  plotAddPairInfo(H_axes,cellInfo);  
end

function addAnnotations(pairUid,pdfFile,xyAreas,chanStr,rowNames,colNames)
    % 
    fontSize = 18;%24
    annotation('textbox',[0.02 0.97 0.05 0.03],'String',pairUid,'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    [~,fn,ext] = fileparts(pdfFile);
    annotation('textbox',[0.15 0.97 0.05 0.03],'String',[fn ext],'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    annotation('textbox',[0.50 0.97 0.05 0.03],'String',xyAreas,'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    h = annotation('textbox',[0.60 0.97 0.10 0.03],'String',chanStr,'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none');
    if contains(chanStr,'Same')
        set(h,'Color','r');
    end
    % conditions / alignNames
    annotation('textbox',[0.01 0.915 0.05 0.05],'String',rowNames{1},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','green')
    annotation('textbox',[0.07 0.895 0.05 0.05],'String',colNames{1},'FontSize',12,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.22 0.895 0.05 0.05],'String',colNames{2},'FontSize',12,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.38 0.895 0.05 0.05],'String',colNames{3},'FontSize',12,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.56 0.895 0.05 0.05],'String',colNames{4},'FontSize',12,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.01 0.485 0.05 0.05],'String',rowNames{2},'FontSize',16,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','red')
end

function doXLabel(H_axis,xLabel)
 fs = 6;
 if ismac
     fs = 8;
 end

 yLim= get(H_axis,'Ylim');
 xPos = mean(get(H_axis,'Xlim')); % center in x
 if strcmp(get(H_axis,'YDir'),'reverse')
    yPos = yLim(2) + range(yLim)*0.05;
 else
    yPos = yLim(1) - range(yLim)*0.05; % 5% below the x-axis
 end
 
 xlabel(xLabel,'Position',[xPos yPos],'VerticalAlignment','top',...
     'HorizontalAlignment','center','FontSize',fs,'FontWeight','bold',...
     'FontAngle', 'italic','Color','black'); 
end

function doYLabel(H_axis,yLabel)
 fs = 6;
 if ismac
     fs = 8;
 end
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
     'HorizontalAlignment','center','FontSize',fs,'FontWeight','bold',...
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
    H_Figure = newFigure();
end

function saveFigAs(fn)

    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',screenposition(3:4),...
        'PaperOrientation','landscape');
    tic
    fprintf('Saving figure to: %s...',fn);
    if ispc
        print(fn,'-dpdf','-painters')
    elseif ismac
        print(fn,'-dpdf','-opengl')
    end
    fprintf('%d\n',toc)
    drawnow
end

function [sortMarkers] = getAlignedSortMarkers(thisTbl)
   % CueOn time is always 3500.
    fx_getFirstSortColIdx = @(tbl) find(strcmp(tbl.Properties.VariableNames,'firstSortByTime'));
    fx_getSecondSortColIdx = @(tbl) find(strcmp(tbl.Properties.VariableNames,'secondSortByTime'));
    alignTimeMinusCueOnTime = thisTbl.alignTime{1} - 3500;
    
    sortMarkers = cell(2,1);
    if ~isempty(fx_getFirstSortColIdx(thisTbl))
        temp = thisTbl{1,fx_getFirstSortColIdx(thisTbl)}{1};
        temp(temp<=0 | isnan(temp)) = -Inf;
        sortMarkers{1} = temp - alignTimeMinusCueOnTime;
    end
    if ~isempty(fx_getSecondSortColIdx(thisTbl))
        temp = thisTbl{1,fx_getSecondSortColIdx(thisTbl)}{1};
        temp(temp<=0 | isnan(temp)) = -Inf;
        sortMarkers{2} =  temp - alignTimeMinusCueOnTime;
    end
end





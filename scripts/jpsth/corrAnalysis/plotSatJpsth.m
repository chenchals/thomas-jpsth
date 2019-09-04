function [] = plotSatJpsth(jpsthPairFile,pdfOutputDir,savePdfFlag)
%PLOTJPSTH Summary of this function goes here


%% Put this in the function....
fx_gSmoothW = @(x,w) smoothdata(x,'gaussian',w,'omitnan');
fx_gSmooth = @(x) fx_gSmoothW(x,20);
conditionPairs = {
    {'Fast','Accurate'};
    {'FastCorrect','AccurateCorrect'};
    {'FastErrorChoice','AccurateErrorChoice'};
    {'FastErrorTiming','AccurateErrorTiming'};
    };
pdfPrefixMap = containers.Map();
pdfPrefixMap(conditionPairs{1}{1}) = 'SAT_';
pdfPrefixMap(conditionPairs{2}{1}) = 'SAT_CORRECT_';
pdfPrefixMap(conditionPairs{3}{1}) = 'SAT_ERROR_CHOICE_';
pdfPrefixMap(conditionPairs{4}{1}) = 'SAT_ERROR_TIMING_';

%%
[~,pdfBaseFile] = fileparts(jpsthPairFile);
for cc = 1:numel(conditionPairs)
    
    conditions = conditionPairs{cc};
    outPdfFile = fullfile(pdfOutputDir, [pdfPrefixMap(conditions{1}) pdfBaseFile '.pdf']);
    jpsthData = load(jpsthPairFile,conditions{:});
    
    %% plot each pair of conditions
    parentFig = getFigHandle();
    H_out = struct();
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);
    offsetsX = [0.005 0.34 0.67]; % for 3 columns
    offsetsY = [0.90 0.45]; % for 2 rows
    startPos = 0.015; % top position of yPsth
    psthH = 0.05; psthW = psthH*3.5;
    gutter = 0.005; % space between plots
    %% compute the min-max for axis scaling

    %#ok<*AGROW>
    for ii = 1:numel(conditions)
        condition = conditions{ii};
        currJpsths = jpsthData.(condition);
        maxSpkPerSecAll(:,ii) = cellfun(@(x) max(fx_gSmooth(x)),[currJpsths.xPsth;currJpsths.yPsth]);
        maxCoinsAll(:,ii) = cellfun(@(x) max(fx_gSmooth(x(:,2))), currJpsths.coincidenceHist);
        minJpsthAll(:,ii) = cellfun(@(x) min(x(:)),currJpsths.normalizedJpsth);
        maxJpsthAll(:,ii) = cellfun(@(x) max(x(:)),currJpsths.normalizedJpsth);
        minXcorrAll(:,ii) =  cellfun(@(x) min(fx_gSmooth(x(:,2))), currJpsths.xCorrHist);
        maxXcorrAll(:,ii) =  cellfun(@(x) max(fx_gSmooth(x(:,2))), currJpsths.xCorrHist);
        minBrodyAll(:,ii) =  cellfun(@(x) min(fx_gSmooth(x)), currJpsths.brodyCovariogram);
        maxBrodyAll(:,ii) =  cellfun(@(x) max(fx_gSmooth(x)), currJpsths.brodyCovariogram);
    end
    maxSpkPerSec = ceil(max(maxSpkPerSecAll(:))*10)/10;
    maxCoins = ceil(max(maxCoinsAll(:))*10)/10;
    minJpsth = floor(min(minJpsthAll(:))*10)/10;
    maxJpsth = ceil(max(maxJpsthAll(:))*10)/10;
    [minXcorr, maxXcorr] = getMinMaxScale(minXcorrAll(:),maxXcorrAll(:));
    [minBrody, maxBrody] = getMinMaxScale(minBrodyAll(:),maxBrodyAll(:));
    
    %% common for all plots
    psthYLims = [0 maxSpkPerSec];
    psthYTicks = [0 maxSpkPerSec];
    psthYTickLabel =  arrayfun(@(x) num2str(x,'%.1f'),psthYTicks','UniformOutput',false);
    
    coinsLims = [-maxCoins maxCoins];
    coinsTicks = [-maxCoins maxCoins];
    coinsTicksLabel =  arrayfun(@(x) num2str(x,'%.1f'),coinsTicks','UniformOutput',false);
    
    normJpsthScale = [minJpsth maxJpsth];
    
    xcorrYLims = [minXcorr maxXcorr];
    xcorrYTicks = [minXcorr maxXcorr];
    xcorrYTickLabel =  arrayfun(@(x) num2str(x,'%.2f'),xcorrYTicks','UniformOutput',false);
    
    brodyYLims = [minBrody maxBrody];
    brodyYTicks = [minBrody maxBrody];
    brodyYTickLabel =  arrayfun(@(x) num2str(x,'%.2f'),brodyYTicks','UniformOutput',false);
    brodySigmaTxtY = range(brodyYLims)/4;
    
    %% plot each condition in a row
    for rowNum = 1:2
        condition = conditions{rowNum};
        axColor = [0.5 0.5 0.5];
        currJpsths = jpsthData.(condition);
        for colNum = 1:3
            alignWinNames = currJpsths.Properties.RowNames;
            % bins for xUnit , yUnit
            psthBins = currJpsths.yPsthBins{colNum};
            psthXLims = [min(psthBins) max(psthBins)];
            psthXTicks = min(psthBins):100:max(psthBins);
            psthXTickLabel =  arrayfun(@(x) num2str(x/1000,'%.1f'),psthXTicks','UniformOutput',false);
            psthXTickLabel(2:end-1) = {''};
            psthXTickLabel(psthXTicks==0) = {'0'};
            % Coincidence Hist
            coinsHist = currJpsths.coincidenceHist{colNum};
            coinsHist(:,1) = psthBins;
            % XCorr hist
            xcorrHist = currJpsths.xCorrHist{colNum};
            xcorrBins = xcorrHist(:,1);
            xcorrXLims = [min(xcorrBins) max(xcorrBins)];
            xcorrXTicks = min(xcorrBins):100:max(xcorrBins);
            xcorrXTickLabel =  arrayfun(@(x) num2str(x/1000,'%.1f'),xcorrXTicks','UniformOutput',false);
            
            % brodyCovariogram (we compute the full lag)
            brodySig = currJpsths.sigma{colNum};
            brodyHist = [xcorrBins currJpsths.brodyCovariogram{colNum}];
            brodyBins = brodyHist(:,1);
            brodyXLims = [min(brodyBins) max(brodyBins)];
            brodyXTicks = min(brodyBins):100:max(brodyBins);
            brodyXTickLabel = arrayfun(@(x) num2str(x/1000,'%.1f'),brodyXTicks','UniformOutput',false);
            brodySigmaTxtX = range(brodyXLims)/4;
            
            %% H_yPsth
            yPsthPos(1) = offsetsX(colNum) + startPos;
            yPsthPos(2) = offsetsY(rowNum) - startPos - psthW;
            yPsthPos(3:4) = [psthH psthW];
            H_out.H_yPsth=axes('parent',parentFig,'position',yPsthPos,'box','on', 'layer','top','Tag','H_yPsth');
            plot(psthBins,fx_gSmooth(currJpsths.yPsth{colNum}));
            view([-90 90])
            annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,90,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0,axColor);
            
            %% H_jpsth
            jpsthPos(1) = yPsthPos(1) + yPsthPos(3) + 3*gutter/aspectRatio;
            jpsthPos(2) = yPsthPos(2);
            jpsthPos(3:4) = [psthW/aspectRatio, psthW];
            %% H_coins1
            coinsPos(1) = jpsthPos(1) + jpsthPos(3) + 3*gutter;
            coinsPos(2) = jpsthPos(2) + jpsthPos(4) - (psthH)*aspectRatio ;
            coinsPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio]; % jpsthPos(3:4);
            H_out.H_coins1=axes('parent',parentFig,'position',coinsPos,'box','on','layer','bottom','Tag','H_coins1');
            area(coinsHist(:,1),fx_gSmooth(coinsHist(:,2)));
            annotateAxis(gca,'y',coinsLims,coinsTicks,coinsTicksLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0,axColor);
            
            % camzoom(sqrt(2));
            % camorbit(-45,0);
            %% H_jpsth
            H_out.H_jpsth=axes('parent',parentFig,'position',jpsthPos,'box','on','layer','top','Tag','H_jpsth');
            imagesc(currJpsths.normalizedJpsth{colNum},normJpsthScale);
            set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
            set(gca,'YDir','normal');
            %% H_colrbar..
            H_out.H_colorbar = colorbar('northoutside');
            colorbarPos = get(H_out.H_colorbar,'Position');
            colorbarPos(2) = jpsthPos(2) + jpsthPos(4) + gutter;
            set(H_out.H_colorbar,'Position',colorbarPos);
            
            %% H_xPsth
            xPsthPos(1) = jpsthPos(1);
            xPsthPos(2) = jpsthPos(2) - 3*gutter - psthH*aspectRatio;
            xPsthPos(3) = psthW/aspectRatio;
            xPsthPos(4) = psthH*aspectRatio;
            H_out.H_xPsth=axes('parent',parentFig,'position',xPsthPos,'box','on','layer','top','Tag','H_xPsth');
            plot(psthBins,fx_gSmooth(currJpsths.xPsth{colNum}));
            annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,0,axColor);
            annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0,axColor);
            
            %% H_xCorr
            xCorrPos(1) = coinsPos(1);
            xCorrPos(2) = coinsPos(2) - (psthH)*aspectRatio - 5*gutter;
            xCorrPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio];
            H_out.H_xCorrHist=axes('parent',parentFig,'position',xCorrPos,'box','on','layer','top','Tag','H_xCorrHist');
            
            area(xcorrHist(:,1),fx_gSmooth(xcorrHist(:,2)));
            annotateAxis(gca,'y',xcorrYLims,xcorrYTicks,xcorrYTickLabel,0,axColor);
            annotateAxis(gca,'x',xcorrXLims,xcorrXTicks,xcorrXTickLabel,0,axColor);
            set(gca,'YAxisLocation','right')
            
            % H_brodyCovariogram (shuffle corrcted cross correlogram)
            xBrodyPos(1) = xCorrPos(1);
            xBrodyPos(2) = xCorrPos(2) - (psthH)*aspectRatio - 5*gutter;
            xBrodyPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio];
            H_out.H_xBrodyHist=axes('parent',parentFig,'position',xBrodyPos,'box','on','layer','top','Tag','H_xBrodyHist');
            
            area(brodyHist(:,1),fx_gSmooth(brodyHist(:,2)));
            hold on
            plot(brodyHist(:,1),1*brodySig,'--c')
            plot(brodyHist(:,1),-1*brodySig,'--c')
            plot(brodyHist(:,1),2*brodySig,'--b')
            plot(brodyHist(:,1),-2*brodySig,'--b')
            plot(brodyHist(:,1),3*brodySig,'--r')
            plot(brodyHist(:,1),-3*brodySig,'--r')
            text(brodySigmaTxtX,brodySigmaTxtY,num2str(max(brodySig),'\\sigma = %.4f'),'Interpreter','tex');
            
            annotateAxis(gca,'y',brodyYLims,brodyYTicks,brodyYTickLabel,0,axColor);
            annotateAxis(gca,'x',brodyXLims,brodyXTicks,brodyXTickLabel,0,axColor);
            set(gca,'YAxisLocation','right')
        end
    end
    drawnow;
    if savePdfFlag
        saveFigAs(outPdfFile);
        delete(parentFig);
    end
end

end

% 1st jpsth plot

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
set(0,'units','pixels');
set(0,'defaulttextfontsize',8,...
    'defaulttextfontname','Arial',...
    'defaultaxesfontsize',8,...
    'defaultaxeslinewidth',0.05);
margin = 10; %pixels
ss=get(0,'ScreenSize');
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

function [minVal, maxVal] = getMinMaxScale(minValues, maxValues)
    result(1) = min(minValues(:));
    result(2) = max(maxValues(:));
    if sum(sign(result)) == 0 % negative and positive
        minVal = -scaleValue(abs(result(1)));
        maxVal = scaleValue(result(2));
    elseif sum(sign(result)) > 0 % only positive
        minVal = 0;
        maxVal = scaleValue(result(2));
    elseif sum(sign(result)) < 0 % only negative
        minVal = -scaleValue(abs(result(1)));
        maxVal = 0;
    end
end

function [tempScale] = scaleValue(temp)
    if temp<0.01
        tempScale = 0.01;
    elseif temp<0.1
        tempScale = 0.1;
    elseif temp<0.5
        tempScale = 0.5;
    else
        tempScale = 1.0;
    end
end



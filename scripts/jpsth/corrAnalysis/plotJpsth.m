%function [] = plotJpsth(jpsthPairFile,varargin)
%PLOTJPSTH Summary of this function goes here

pairDir = 'dataProcessed/analysis/JPSTH/jpsth_SEF-FEF';
jpsthPairFile = fullfile(pairDir,'JPSTH-PAIR_0027.mat');
jpsthData = load(jpsthPairFile);
%%
conditions = {'FastErrorTiming','AccurateErrorTiming'};
fx_gSmoothW = @(x,w) smoothdata(x,'gaussian',w,'omitnan');
fx_gSmooth = @(x) fx_gSmoothW(x,20);
parentFig = getFigHandle();
H_out = struct();
ss = get(0,'ScreenSize');
aspectRatio = ss(3)/ss(4);
maxRows = 2; % max rows of jpsths
maxCols = 3;% max jpsths in a row
offsetsX = [0.005 0.34 0.67]; % for 3 columns
offsetsY = [0.94 0.48]; % for 2 rows
startPos = 0.015; % top position of yPsth
psthH = 0.05; psthW = psthH*3.5;
gutter = 0.005; % space between plots
%% compute the min-max for axis scaling
%#ok<*SAGROW>
for ii = 1:numel(conditions)
    condition = conditions{ii};
    currJpsths = jpsthData.(condition);
    maxSpkPerSecAll(:,ii) = cellfun(@(x) max(fx_gSmooth(x)),[currJpsths.xPsth;currJpsths.yPsth]);
    maxCoinsAll(:,ii) = cellfun(@(x) max(fx_gSmooth(x(:,2))), currJpsths.coincidenceHist);
    minJpsthAll(:,ii) = cellfun(@(x) min(x(:)),currJpsths.normalizedJpsth);
    maxJpsthAll(:,ii) = cellfun(@(x) max(x(:)),currJpsths.normalizedJpsth);
    minXcorrAll(:,ii) =  cellfun(@(x) min(fx_gSmooth(x(:,2))), currJpsths.xCorrHist);
    maxXcorrAll(:,ii) =  cellfun(@(x) max(fx_gSmooth(x(:,2))), currJpsths.xCorrHist);
end
maxSpkPerSec = ceil(max(maxSpkPerSecAll(:))*10)/10;
maxCoins = ceil(max(maxCoinsAll(:))*10)/10;
minJpsth = floor(min(minJpsthAll(:))*10)/10;
maxJpsth = ceil(max(maxJpsthAll(:))*10)/10;

temp = abs(min(minXcorrAll(:)));
if temp<0.01
    minXcorr = -0.01;
elseif temp<0.1
    minXcorr = -0.1;
elseif temp<0.5
    minXcorr = -0.5;
else
    minXcorr = -1.0;
end

temp = max(maxXcorrAll(:));
if temp<0.01
    maxXcorr = 0.01;
elseif temp<0.1
    maxXcorr = 0.1;
elseif temp<0.5
    maxXcorr = 0.5;
else
    maxXcorr = 1.0;
end

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


%% plot each condition in a row
for rowNum = 1:2
    condition = conditions{rowNum};
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
        
        %% H_yPsth
        yPsthPos(1) = offsetsX(colNum) + startPos;
        yPsthPos(2) = offsetsY(rowNum) - startPos - psthW;
        yPsthPos(3:4) = [psthH psthW];
        H_out.H_yPsth=axes('parent',parentFig,'position',yPsthPos,'box','on','layer','top','Tag','H_yPsth');
        plot(psthBins,fx_gSmooth(currJpsths.yPsth{colNum}));
        view([-90 90])     
        annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,90);
        annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0);

        %% H_jpsth
        jpsthPos(1) = yPsthPos(1) + yPsthPos(3) + 2*gutter;
        jpsthPos(2) = yPsthPos(2);
        jpsthPos(3:4) = [psthW/aspectRatio, psthW];
        %% H_coins1
        coinsPos(1) = jpsthPos(1) + jpsthPos(3) + 3*gutter;
        coinsPos(2) = jpsthPos(2) + jpsthPos(4) - (psthH + psthH/4)*aspectRatio ;
        coinsPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio]; % jpsthPos(3:4);       
        H_out.H_coins1=axes('parent',parentFig,'position',coinsPos,'box','on','layer','bottom','Tag','H_coins1');
        area(coinsHist(:,1),fx_gSmooth(coinsHist(:,2)));        
        annotateAxis(gca,'y',coinsLims,coinsTicks,coinsTicksLabel,0);
        annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0);

        % camzoom(sqrt(2));
        % camorbit(-45,0);
        %% H_jpsth
        H_out.H_jpsth=axes('parent',parentFig,'position',jpsthPos,'box','on','layer','top','Tag','H_jpsth');
        imagesc(currJpsths.normalizedJpsth{colNum},normJpsthScale);
        set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
        set(gca,'YDir','normal');
        %% H_colrbar..
        
        %% H_xPsth
        xPsthPos(1) = jpsthPos(1);
        xPsthPos(2) = jpsthPos(2) - 2*gutter - psthH*aspectRatio;
        xPsthPos(3) = psthW/aspectRatio;
        xPsthPos(4) = psthH*aspectRatio;
        H_out.H_xPsth=axes('parent',parentFig,'position',xPsthPos,'box','on','layer','top','Tag','H_xPsth');
        plot(psthBins,fx_gSmooth(currJpsths.xPsth{colNum}));     
        annotateAxis(gca,'y',psthYLims,psthYTicks,psthYTickLabel,0);
        annotateAxis(gca,'x',psthXLims,psthXTicks,psthXTickLabel,0);
        
        %% H_xCorr
        xCorrPos(1) = coinsPos(1);
        xCorrPos(2) = coinsPos(2) - (psthH + psthH/4)*aspectRatio - 5*gutter;
        xCorrPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio];
        H_out.H_xCorrHist=axes('parent',parentFig,'position',xCorrPos,'box','on','layer','top','Tag','H_xPsth');

        area(xcorrHist(:,1),fx_gSmooth(xcorrHist(:,2)));        
        annotateAxis(gca,'y',xcorrYLims,xcorrYTicks,xcorrYTickLabel,0);
        annotateAxis(gca,'x',xcorrXLims,xcorrXTicks,xcorrXTickLabel,0);  
        set(gca,'YAxisLocation','right')
        
    end
end



% 1st jpsth plot

function annotateAxis(H_axis,xOry,lims_,ticks_,labels_,rotDeg)
    if xOry=='y'
        set(H_axis,'ylim',lims_,'ytick',ticks_,'YTickLabel',labels_,...
            'YTickLabelRotation',rotDeg,'ycolor',[0.5 0.5 0.5]);
    elseif xOry=='x'
        set(H_axis,'xlim',lims_,'xtick',ticks_,'XTickLabel',labels_,...
            'XTickLabelRotation',rotDeg,'xcolor',[0.5 0.5 0.5]);
        set(H_axis,'XGrid','on','GridLineStyle','--', 'GridColor', [0.3 0.3 0.3])
        line([0 0],lims_,'Color',[0.1 0.1 0.1],'LineStyle',':')
    end
end

%%
function [H_Figure] = getFigHandle()
    scale = 0.7;
    set(0,'units','pixels');
    set(0,'defaulttextfontsize',5*scale,...
        'defaulttextfontname','Arial',...
        'defaultaxesfontsize',5*scale,...
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




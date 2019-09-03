%function [] = plotJpsth(jpsthPairFile,varargin)
%PLOTJPSTH Summary of this function goes here

% pairDir = 'dataProcessed/analysis/JPSTH/jpsth_SEF-FEF';
% jpsthPairFile = fullfile(pairDir,'JPSH-PAIR_0008.mat');
% pairData = load(jpsthPairFile);
%%
conditions = {'Fast','Accurate'};

H_Figure = getFigHandle();
H_out = struct();
ss = get(0,'ScreenSize');
aspectRatio = ss(3)/ss(4);
maxRows = 2; % max rows of jpsths
maxCols = 3;% max jpsths in a row
offsetsX = [0 0.32 0.65]; % for 3 columns
offsetsY = [0.95 0.48]; % for 2 rows
startPos = [0.015]; % top position of yPsth
psthH = 0.05; psthW = psthH*4;
gutter = 0.002; % space between plots

for rowNum = 1:2
    condition = conditions{row};
    for colNum = 1:3
        %        H_Plots{row,col} =  getJsthAxesHandles(H_Figure,row,col);
        
        % H_yPsth
        yPsthPos(1) = offsetsX(colNum) + startPos;
        yPsthPos(2) = offsetsY(rowNum) - startPos - psthW;
        yPsthPos(3:4) = [psthH psthW];
        H_out.H_yPsth=axes('parent',parentFig,'position',yPsthPos,'box','on','layer','top','Tag','H_yPsth');
        plot(1:100,rand(1,100));
        view([-90 90])
        set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
        % H_jpsth
        jpsthPos(1) = yPsthPos(1) + yPsthPos(3) + gutter;
        jpsthPos(2) = yPsthPos(2);
        jpsthPos(3:4) = [psthW/aspectRatio, psthW];
        % H_coins1
        coinsPos(1) = jpsthPos(1) + jpsthPos(3) + gutter;
        coinsPos(2) = jpsthPos(2) + jpsthPos(4) - (psthH + psthH/4)*aspectRatio ;
        coinsPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio]; % jpsthPos(3:4);
        H_out.H_coins1=axes('parent',parentFig,'position',coinsPos,'box','on','layer','bottom','Tag','H_coins1');
        area(1:100,rand(1,100)-0.5)
        set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
        ylim([-1 1])
        % camzoom(sqrt(2));
        % camorbit(-45,0);
        % H_jpsth
        H_out.H_jpsth=axes('parent',parentFig,'position',jpsthPos,'box','on','layer','top','Tag','H_jpsth');
        imagesc(rand(500,500))
        set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
        set(gca,'YDir','normal');
        % H_colrbar..
        % H_xPsth
        xPsthPos(1) = jpsthPos(1);
        xPsthPos(2) = jpsthPos(2) - gutter - psthH*aspectRatio;
        xPsthPos(3) = psthW/aspectRatio;
        xPsthPos(4) = psthH*aspectRatio;
        H_out.H_xPsth=axes('parent',parentFig,'position',xPsthPos,'box','on','layer','top','Tag','H_xPsth');
        plot(1:100,rand(1,100));
        set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
        % H_xCorr
        xCorrPos(1) = coinsPos(1);
        xCorrPos(2) = coinsPos(2) - (psthH + psthH/4)*aspectRatio - 5*gutter;
        xCorrPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio];
        H_out.H_xCorrHist=axes('parent',parentFig,'position',xCorrPos,'box','on','layer','top','Tag','H_xPsth');
        plot(1:100,rand(1,100));
        set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
        
        
        
        
        
    end
end



% 1st jpsth plot

%%
function [H_Figure] = getFigHandle()
    scale = 0.7;
    set(0,'units','pixels');
    set(0,'defaulttextfontsize',6*scale,...
        'defaulttextfontname','Arial',...
        'defaultaxesfontsize',6*scale,...
        'defaultaxeslinewidth',0.05);
    margin = 20; %pixels
    ss=get(0,'ScreenSize');
    FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
    %Main figure window
    H_Figure=figure('Position',FigPos,...
        'color',[1 1 1],'numbertitle','off','renderer','painters',...
        'renderermode','manual','menubar','none','name','',...
        'Tag','H_Figure');
    orient landscape
    set(H_Figure,'units','normalized')
end

function [H_out] = getJsthAxesHandles(parentFig, rowNum, colNum, jpsthTable)
  H_out = struct();
  ss = get(0,'ScreenSize');
  aspectRatio = ss(3)/ss(4);
  maxRows = 2; % max rows of jpsths
  maxCols = 3;% max jpsths in a row
  offsetsX = [0 0.32 0.65]; % for 3 columns
  offsetsY = [0.95 0.48]; % for 2 rows
  startPos = [0.015]; % top position of yPsth
  psthH = 0.05; psthW = psthH*4;
  gutter = 0.002; % space between plots
  % H_yPsth
  yPsthPos(1) = offsetsX(colNum) + startPos;
  yPsthPos(2) = offsetsY(rowNum) - startPos - psthW;
  yPsthPos(3:4) = [psthH psthW]; 
  H_out.H_yPsth=axes('parent',parentFig,'position',yPsthPos,'box','on','layer','top','Tag','H_yPsth');
   plot(1:100,rand(1,100));
   view([-90 90])
   set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
 % H_jpsth
 jpsthPos(1) = yPsthPos(1) + yPsthPos(3) + gutter;
 jpsthPos(2) = yPsthPos(2);
 jpsthPos(3:4) = [psthW/aspectRatio, psthW];
 % H_coins1
 coinsPos(1) = jpsthPos(1) + jpsthPos(3) + gutter;
 coinsPos(2) = jpsthPos(2) + jpsthPos(4) - (psthH + psthH/4)*aspectRatio ;
 coinsPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio]; % jpsthPos(3:4);
 H_out.H_coins1=axes('parent',parentFig,'position',coinsPos,'box','on','layer','bottom','Tag','H_coins1');
 area(1:100,rand(1,100)-0.5)
 set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
 ylim([-1 1])
 % camzoom(sqrt(2));
 % camorbit(-45,0);
 % H_jpsth
 H_out.H_jpsth=axes('parent',parentFig,'position',jpsthPos,'box','on','layer','top','Tag','H_jpsth');
 imagesc(rand(500,500))
 set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
 set(gca,'YDir','normal');
 % H_colrbar..
 % H_xPsth
 xPsthPos(1) = jpsthPos(1);
 xPsthPos(2) = jpsthPos(2) - gutter - psthH*aspectRatio;
 xPsthPos(3) = psthW/aspectRatio;
 xPsthPos(4) = psthH*aspectRatio;
 H_out.H_xPsth=axes('parent',parentFig,'position',xPsthPos,'box','on','layer','top','Tag','H_xPsth');
 plot(1:100,rand(1,100));
 set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
 % H_xCorr
 xCorrPos(1) = coinsPos(1);
 xCorrPos(2) = coinsPos(2) - (psthH + psthH/4)*aspectRatio - 5*gutter;
 xCorrPos(3:4) = [psthW/aspectRatio, psthH*aspectRatio]; 
 H_out.H_xCorrHist=axes('parent',parentFig,'position',xCorrPos,'box','on','layer','top','Tag','H_xPsth');
 plot(1:100,rand(1,100));
 set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);


end

function [rotX,rotY] = getRotatedVals(xVals,yVals)
% From: https://stackoverflow.com/questions/29780903/2d-body-transformation-and-rotation-in-matlab
zVals = ones(1,numel(xVals));
xVals = xVals(:)';
yVals = yVals(:)';
transOrig = @(x,y,z) repmat([x;y;z;],[1,numel(x)]);
rotOrig =  @(x, y, theta) [
    cosd(theta), -sind(theta), x;
    sind(theta), cosd(theta), y;
    0,            0,           1];
rotMat1 = rotOrig(0,0,45) ...
    * ([xVals;yVals;zVals] - transOrig(floor(numel(xVals)/2),0,0) )...
    + transOrig(floor(numel(xVals)/2),0,0);
rotMat = rotOrig(0,0,45) * ([xVals;yVals;zVals])
rotX = rotMat(1,:)';
rotY = rotMat(2,:)';

end

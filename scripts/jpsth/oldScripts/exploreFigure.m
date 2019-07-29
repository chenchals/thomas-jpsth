% draw jpsth fig




% input is spiketimes or binned spike times
%figDir = '/Volumes/schalllab/Users/Chenchal/JPSTH/FigsBinWidth_5/PAIR_0226_E20130829001-RH_SEARCH_DSP16a_SEF_DSP17a_SC.mat';%
figDir = '/Volumes/schalllab/Users/Chenchal/JPSTH/FigsBinWidth_5';%
files = dir(fullfile(figDir,'PAIR*.mat'));
fns = {files.name}';
files = strcat({files.folder}',filesep,{files.name}');
conditions = whos('-file',files{1});
conditions = {conditions.name}';
conditions(strcmpi(conditions,'cellPairInfo'))=[];
alignedEvents = {'CueOn';'SaccadePrimary';'RewardOn'};
jpMinMax = struct();
fx_minmax = @(x) [min(x(:)) max(x(:))];
for ii = 1:numel(files)
    currFile = files{ii};
    [~,pairName] = fileparts(currFile);
    jpMinMax(ii).pairName = pairName;
    tic
    pairConds = load(currFile);
    cellPairInfo = pairConds.cellPairInfo;
    for cc = 1:numel(conditions)
        condition = conditions{cc};
        pair = pairConds.(condition);
        if isempty(pair)
            jpMinMax(ii).(condition) = nan(3,2);
        else

            jpsth = pair.normalizedJpsth;
            jpMinMax(ii).(condition) = cell2mat(cellfun(fx_minmax,jpsth,'UniformOutput',false));
            jpMinMax(ii).(['grand_' condition])= fx_minmax(jpMinMax(ii).(condition));

            %imagesc(flipud(jpsth));
            
        end
        title(pairName,'Interpreter','none')
    end % each condition
    toc
end

%% Plot coincidence
  fontscale=1.0;
  linescale=0.05;
  coins = Z.AccurateCorrect.coincidenceHist{1};
  binWidth = Z.AccurateCorrect.coincidenceBins(1);
  coins(:,1) = coins(:,1)*binWidth;
  xlims = minmax(coins(:,1)');
  ylims =[-max(abs(coins(:,2))) max(abs(coins(:,2)))];

   set(bar(coins(:,1),coins(:,2),1),'facecolor',[1.0 0.0 1.0],'edgecolor','none')
   set(gca,'xlim',xlims)
   set(gca,'ylim',ylims)
   set(gca,'xticklabel',[],'yticklabel',[])
   set(line(xlims,[0 0]),'color',[0 0 0],'linewidth',linescale*1.2);
   set(line([max(xlims) max(xlims)],[0 max(ylims)]),'color',[0 0 0],'linewidth',linescale*1.0);
   set(line([max(xlims) max(xlims)-(4*binWidth)],[max(ylims) max(ylims)]),'color',[0 0 0],'linewidth',linescale*1.0);
   
   text(max(xlims)*0.96,max(ylims)*0.96,num2str(max(ylims),'%2.3f '),...
      'horizontalalignment','right','verticalalignment','bottom',...
      'fontname','Arial','fontsize',fontscale*6,'fontweight','normal','rotation',45);
   pbaspect([3 1 1])
   camzoom(sqrt(2));
   camorbit(-45,0)
   set(gca,'layer','bottom')
   set(gca,'Box','off')
   
   %% Plot xCorr..
   
   correl = Z.AccurateCorrect.xCorrHist{1};
   xVals = correl(:,1);
   yVals = correl(:,2);
   binWidth = 5;
   AboveMain=5*binWidth; %round(XFiresBeforeY/BINSIZE);
   BelowMain=-5*binWidth; %round(XFiresAfterY/BINSIZE);

   
      %always have a total of 9 ticks for x axis (including 0)
   xc_tick=[];step=[];
   step=(range(xVals))/8;
   xc_tick=([-4:1:4])'.*step;
   temp=yVals;
   minmax=[min(temp) max(temp)];
   y_axis=[];
   y_axis=minmax;
   if max(temp)<0 ,y_axis(2)=0;end
   if min(temp)>0 ,y_axis(1)=0;end
   yscale=0;
   if abs(min(temp))> abs(max(temp))
      yscale=min(temp);
   else
      yscale=max(temp);
   end
      
   set(bar(correl(:,1),correl(:,2),1),'facecolor',[1 0.25 1],'edgecolor','none')
   
   set(gca,'yticklabel',[],'xticklabel',[],'xtick',xc_tick)
   set(gca,'xlim',[-range(xVals) range(xVals)],...
      'ylim',y_axis);
   set(line([-range(xVals) range(xVals)],[0 0]),'color',[0 0 0],'linewidth',linescale*1.0);
   set(line([-range(xVals) -range(xVals)],[0 yscale]),'color',[0 0 0],'linewidth',linescale*1.0);
   set(line([-range(xVals) -range(xVals)+4*binWidth],[yscale yscale]),'color',[0 0 0],'linewidth',linescale*1.0);
   text(-range(xVals)*0.9,yscale,num2str(yscale,'%2.3f '),...
      'horizontalalignment','left','verticalalignment','middle',...
      'fontname','Arial','fontsize',fontscale*6,'fontweight','normal','rotation',-45);
   set(line([(AboveMain+0.5)*binWidth (BelowMain+0.5)*binWidth], [yscale*1.05 yscale*1.05] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([(AboveMain+0.5)*binWidth (AboveMain+0.5)*binWidth], [yscale*1.05 yscale*0.99] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([(BelowMain+0.5)*binWidth (BelowMain+0.5)*binWidth], [yscale*1.05 yscale*0.99] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([0 0], [max(yVals)*1.07 max(yVals)*1.12] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([0 0], [0 yscale*1.15]),'color',[0 0 0],'linewidth',linescale*1.0)
   pbaspect([5 1 1]);
   camzoom(sqrt(2))
   camorbit(45,0)
   set(gca,'Box','off')

   
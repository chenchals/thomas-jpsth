%% For spike distribution plot
figure
subplot(3,1,1)
plot(-600:100,spikeCorr.rho_pval_200ms{1}(:,1),'LineWidth',1.5,'LineStyle','--','Color','b')
hold on
plot(-200:400,spikeCorr.rho_pval_200ms{2}(:,1),'LineWidth',1.5,'LineStyle',':','Color','r')
legend({'signal corr: Baseline','signal Corr: Visual'},'Location','northwest')
xlim([-600 400])
title('Moving window of 200ms - Raw counts' )
ylabel('\rho','Interpreter','tex')
xlabel('Aligned on ''Cue on''')

subplot(3,1,2)
plot(-600:100,spikeCorr.rho_pval_200ms_Z_baseline{1}(:,1),'LineWidth',1.5,'LineStyle','--','Color','b')
hold on
plot(-200:400,spikeCorr.rho_pval_200ms_Z_baseline{2}(:,1),'LineWidth',1.5,'LineStyle',':','Color','r')
xlim([-600 400])
legend({'noise corr: Baseline','noise corr: Visual'},'Location','northwest')
title('Moving window of 200ms - Z-scored using Baseline mean and sigma' )
ylabel('\rho','Interpreter','tex')
xlabel('Aligned on ''Cue on''')

subplot(3,1,3)
plot(-600:100,spikeCorr.rho_pval_200ms_Z_trial{1}(:,1),'LineWidth',1.5,'LineStyle','--','Color','b')
hold on
plot(-200:400,spikeCorr.rho_pval_200ms_Z_trial{2}(:,1),'LineWidth',1.5,'LineStyle',':','Color','r')
xlim([-600 400])
legend({'noise corr: Baseline','noise corr: Visual'},'Location','northwest')
title('Moving window of 200ms - Z-scored using Trial mean and sigma' )
ylabel('\rho','Interpreter','tex')
xlabel('Aligned on ''Cue on''')

%%
figure

binStr = {'10ms','50ms','100ms','200ms'};
lineColors = {'b','g','c','r'};
epochNameIdx = 2; % 1-4 = baseline, visual, postSacc, postRew
epochNames = {'Baseline','Visual','PostSaccade','PostReward'}; % 1= baseline
alignedOnNames = {'''Cue On''', '''Cue On''', '''Saccade Primary''', '''Reward Time''', };
xVals = {-600:100,-200:400,-100:500,-200:700}; % 1=[-600:100],2=[-200:400]
for rz = 1:3   
    subplot(3,1,rz)
    hold on
    grid on
    x = xVals{epochNameIdx};
    colNameSuffix ='';
    epochName = epochNames{epochNameIdx};
    
    ylabel('\rho','Interpreter','tex')
    xlabel(['Aligned on ' alignedOnNames{epochNameIdx}])
    title(['\fontsize{15} ' epochName ' \fontsize{12} Signal correlation - raw counts'],'Interpreter','tex');
    if rz==2
        colNameSuffix ='_Z_baseline';
        title(['\fontsize{15} ' epochName ' \fontsize{12} Noise correlation - Z-scored with Baseline mean and sigma'],'Interpreter','tex');
    elseif rz==3
        colNameSuffix ='_Z_trial';
        title(['\fontsize{15} ' epochName ' \fontsize{12} Noise correlation - Z-scored with Trial mean and sigma'],'Interpreter','tex');
    end
    for b = 1:4
        colName = ['rho_pval_' binStr{b} colNameSuffix];
        rho = spikeCorr.(colName){epochNameIdx}(:,1);
        pval = spikeCorr.(colName){epochNameIdx}(:,2);
        plot(x,rho,'LineWidth',0.05,'LineStyle','-','Color',lineColors{b})
        rho(pval>0.05) = NaN;
        plot(x,rho,'LineWidth',3.5,'LineStyle','-','Color',lineColors{b},'HandleVisibility','off')        
    end
    legend(strcat(binStr,' moving-win.'),'Location','northwest')
end

%%
figure
subplot(2,1,1)
plot(spikeCorr.xBaselineMeanStd{1}(:,1))
hold on
plot(spikeCorr.xTrialMeanStd{1}(:,1))

subplot(2,1,2)
plot(spikeCorr.yBaselineMeanStd{1}(:,1))
hold on
plot(spikeCorr.yTrialMeanStd{1}(:,1))


xVar='xSpkCount_win';
yVar='ySpkCount_win';
rpVar='rho_pval_static';

X=spkCorr.(xVar){1};
Y=spkCorr.(yVar){1};
rho_pval=spkCorr.(rpVar){1};

fx_scale = @(x) (x-min(x))./range(x);
Y = fx_scale(Y);
X = fx_scale(X);

 [b,bint] = regress(Y,X);
 xval = min(X):0.01:max(X);
 yhat = b(1)*xval;
 ylow = bint(1)*xval;
 yupp = bint(2)*xval;
 figure
 scatter(X,Y,'*');
 hold on;
 plot(xval,ylow,'b-.');
 plot(xval,yupp,'r-.');
 plot(xval,yhat,'b','linewidth',2);
 
 
 %% histogram...
 figure
 subplot(1,2,1)
 rawDat = double(currDat.rhoRaw);
 minMaxRho = round(minmax(rawDat(:)'),1); % to 1st decimal
 binWidth = 0.05;
 bins = minMaxRho(1):binWidth:minMaxRho(2);
 y = histcounts(rawDat,bins);
 y(y==0) = NaN;
 ypercent = y.*100/nansum(y);
 x = bins(1:end-1) + binWidth/2;
 h = bar(x,ypercent,'BarWidth',1);
 hold on
 h1 = text(double(x),ypercent,num2str(y(:),'%d'),...
     'VerticalAlignment','bottom','HorizontalAlignment','center',...
     'FontSize',10); 
 ylim([0 max(ypercent)+2])
 ylabel('Percentage of Pairs')
 xlabel('Spike Count Correlation (r_{sc})','Interpreter','tex')
 yrange = range(get(gca,'ylim'));
 xrange = range(get(gca,'xlim'));
 % annotate
 scatter(nanmean(rawDat),yrange*0.9,100,'v','filled','MarkerEdgeColor','r','MarkerFaceColor','r'); 
 line([nanmean(rawDat),nanmean(rawDat)],[0 yrange*0.9],'color','r','LineStyle','--','LineWidth',1);
 text(nanmean(rawDat),yrange*0.95,sprintf('%0.2f (%d)',nanmean(rawDat),numel(rawDat)),'FontSize',12,'FontWeight','bold');
 hold off
 subplot(1,2,2)
 rawxDat = double(currDat.XY_Dist);
 rawyDat = double(currDat.rhoRaw);
 minMaxDist = round(minmax(rawxDat(:)'),1); % to 1st decimal
 binWidth = 1.0;%mm
 %bins = (-binWidth:binWidth:minMaxDist(2))'; 
 bins = (minMaxDist(1):binWidth:ceil(minMaxDist(2)))';
 y = arrayfun(@(x) mean(rawyDat(rawxDat>=bins(x) & rawxDat<bins(x+1))), (1:numel(bins)-1)');
 ystd = arrayfun(@(x) std(rawyDat(rawxDat>=bins(x) & rawxDat<bins(x+1))), (1:numel(bins)-1)');
 ycounts = arrayfun(@(x) sum(rawxDat>=bins(x) & rawxDat<bins(x+1)), (1:numel(bins)-1)');
 x = bins(1:end-1);
 h = bar(x,y,'BarWidth',1);
 hold on
he = errorbar(x,y,ystd,'Marker','o','MarkerSize',10,'LineStyle',...
    'none','LineWidth',1.5,'MarkerFaceColor','w','MarkerEdgeColor','k') 
%set(he,'YPositiveDelta',[])
set(gca,'YLim',[0 max(y)*1.1])
  %mustdn = [num2str(y(:),'%0.2f\\pm') num2str(ystd(:),'%0.2f') repmat(' ',numel(y),1)  num2str(ycounts(:),'n=%d')]
  mustdn = [num2str(ycounts(:),'%d')]
      h1 = text(x+0.1,y+0.005,mustdn,'VerticalAlignment','top','HorizontalAlignment','left',...
         'FontSize',10); 

 
%% grouped error bar plot

a=[0,1,0,0;
4,3,2,1;
2,2,1,3;
1,0,0,0];
b=[0,1,0,0;
1,2,1,1;
1,1,1,2;
1,0,0,0];
ctrs = 1:4;
data = a;
figure(1)
hBar = bar(ctrs, data);
for k1 = 1:size(a,1)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, b, '.r')
hold off


%% cheat for annotation
xOff = 0.01;
yOff = [0.8 0.6 0.3];
h = axes('Parent',gcf,'Position',[xOff yOff(2) 0.005 0.005],'Box','off','Visible','off')


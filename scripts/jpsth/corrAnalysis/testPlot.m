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


figure
subplot(2,1,1)
plot(spikeCorr.xBaselineMeanStd{1}(:,1))
hold on
plot(spikeCorr.xTrialMeanStd{1}(:,1))

subplot(2,1,2)
plot(spikeCorr.yBaselineMeanStd{1}(:,1))
hold on
plot(spikeCorr.yTrialMeanStd{1}(:,1))




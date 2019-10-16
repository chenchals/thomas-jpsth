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
        epoch = alignedNames{an};        
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
        end
    end
    
end


function [] = plotRscScatterStats(x,y,xmu,ymu,ystd,sig05,sig01)
    gca
    col_k = [0 0 0]; col_m = [1 0 1]; col_c = [0 1 1]; col_b=[0 0 1]; col_r = [1 0 0];
    scatter(x,y,30,'o','filled','MarkerFaceColor',col_k,'MarkerFaceAlpha',0.2)
    hold on
    xlim([-0.5 10])
    ylim([-1 1])
    scatter(x(sig05),y(sig05),30,'o','filled','MarkerFaceColor',col_c,'MarkerFaceAlpha',0.5)
    scatter(x(sig01),y(sig01),30,'o','filled','MarkerFaceColor',col_m,'MarkerFaceAlpha',0.5)
    % draw the mean line
    h = plot(xmu,ymu,'Color',col_b,'Marker','s','MarkerFaceColor',col_b, 'LineWidth',1);
    h_var = fill([xmu fliplr(xmu)], [ymu+ystd fliplr(ymu-ystd)],col_b,'FaceAlpha',0.1,'LineStyle','none');
    set([h_var,h],'HandleVisibility','off');
    hold off
end

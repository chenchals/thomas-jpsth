% Compare cross area and within area spike count correlations
% Use static windows
fn = 'dataProcessed/analysis/spkCorr/spkCorrAllPairsStatic.mat';
spkCorr = load(fn);
outFn = 'dataProcessed/analysis/spkCorr/Summary_SpkCorr_AcrossAreas.mat';

% group stats for the following areas and align windows
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
% Prune to columns of interest for each pairArea
useCols = {
    'pairAreas'
    'nTrials'
    'XY_Dist'
    'condition'
    'alignedName'
    'alignedEvent'
    'rho_pval_win'
    'rhoRaw'
    'pvalRaw'
    'rhoZTrial'
    'pvalZTrial'
    'signifRaw_05'
    'signifRaw_01'
    'signifZTrial_05'
    'signifZTrial_01'
    };
compareConditions = {
    'AccurateCorrect' 'FastCorrect'
    'AccurateErrorChoice' 'FastErrorChoice'
    'AccurateErrorTiming' 'FastErrorTiming'
    };
% aggregate required fields into a single table
spkCorrAll = table();
for pa =1:numel(pairAreas)
    pairArea = pairAreas{pa};
    spkCorrAll = [spkCorrAll; spkCorr.(pairArea)(:,useCols)];
end
% use absolute values
spkCorrAll.rhoRaw = abs(spkCorrAll.rhoRaw);
spkCorrAll.rhoZTrial = abs(spkCorrAll.rhoZTrial);
spkCorrAll.baseCondition = regexprep(spkCorrAll.condition,'(Accurate)|(Fast)','');
spkCorrAll.trialType = regexprep(spkCorrAll.condition,'(Correct)|(Error.*)','');

% get stats
spkCorrAllStats = grpstats(spkCorrAll(:,{'baseCondition','trialType','pairAreas','condition','alignedName','rhoRaw','rhoZTrial'}),...
        {'condition','baseCondition','trialType','alignedName','pairAreas'},{'mean','std','min','max'});


%% Plot the stats
baseConditions = {'Correct','ErrorChoice','ErrorTiming'};
alignedNames = {'Baseline','Visual','PostSaccade','PostReward'};

H_plots = getPlotHandles();
H_plt_Idx = 0;
for co = 1:numel(alignedNames)
    alignedName = alignedNames{co};
    for ro = 1:numel(baseConditions)
        baseCond = baseConditions{ro};
        tblRows = strcmp(spkCorrAllStats.alignedName,alignedName) & ...
                  strcmp(spkCorrAllStats.baseCondition,baseCond);
        datAreas = reshape(spkCorrAllStats.pairAreas(tblRows),numel(pairAreas),[]);    
        datMean = reshape(spkCorrAllStats.mean_rhoRaw(tblRows),numel(pairAreas),[]);
        datStd = reshape(spkCorrAllStats.std_rhoRaw(tblRows),numel(pairAreas),[]);        
        % plot it
        H_plt_Idx = H_plt_Idx + 1;
        plotGroupedBarsWithErrBars(H_plots(H_plt_Idx),datAreas,datMean,datStd);        
    end
end

addAnnotations(outFn,baseConditions,alignedNames);
  
saveFigPdf(outFn);
    
    
%% Other functions  

function addAnnotations(pdfFile,rowNames,colNames)
    % 
    fontSize = 18;%24
    [~,fn,ext] = fileparts(pdfFile);
    annotation('textbox',[0.15 0.97 0.05 0.03],'String',[fn ext],'FontSize',fontSize,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','Interpreter','none')
    % conditions / alignNames
    annotation('textbox',[0.14 0.92 0.05 0.05],'String',colNames{1},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.38 0.92 0.05 0.05],'String',colNames{2},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.62 0.92 0.05 0.05],'String',colNames{3},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    annotation('textbox',[0.86 0.92 0.05 0.05],'String',colNames{4},'FontSize',14,'FontWeight','bold','FontAngle','italic','FitBoxToText','on','EdgeColor','none','color','black')
    % row names    
    xOff = 0.02;
    yOff = [0.8 0.5 0.2];
    for ii = 1:numel(rowNames)
         axes('Parent',gcf,'Position',[xOff yOff(ii) 0.005 0.005],'Box','off','Visible','off')
         text(0,0,rowNames{ii},'FontSize',14,'FontWeight','bold','Rotation',90)
    end
end


function [] = plotGroupedBarsWithErrBars(H_axes,xData,yMean,yStd)
    axes(H_axes);
    %x = categorical(xData(:,1));
    x = 1:3;
    hBar = bar(x,yMean,'FaceAlpha',0.4);
    set(hBar(1),'FaceColor',[1 0 0]); set(hBar(2),'FaceColor',[0 1 0]);
    set(gca,'xticklabels',xData(:,1),'TickLabelInterpreter','none');
    ylabel('Spike count corr.') 
    hold on
    for k1 = 1:size(yMean,2)
        ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    end
    he = errorbar(ctr, yMean', yStd','Marker','o','MarkerSize',10,'LineStyle',...
          'none','LineWidth',1.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');
    hold off

end


function [H_plots] = getPlotHandles()
    %% for spike corr comparision among areas
    % 3 rows : correct, errorChoice, errorTiming
    % 4 cols : baseline, visual, postSaccade, postReward
    ss = [20 20 1500 990];
    margin = 10;
    FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
    %Main figure window
    H_Figure=figure('Position',FigPos,...
        'color',[1 1 1],'numbertitle','off','renderer','painters',...
        'renderermode','manual','menubar','none',...
        'Tag','H_Figure');
    orient landscape
    set(H_Figure,'units','normalized')
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);

    nRows = 3;
    nCols = 4;

    pltH = 0.25;
    pltW = 0.20;
    gutter = 0.04;
    offsetsX = 0.06:(pltW+gutter):1-gutter; % for 4 column-starts
    offsetsY = 0.95-pltH:-(pltH+gutter+0.01):gutter; % for 3 row-starts

    pltCount = 0;
    for cols = 1:nCols
        pos(1) = offsetsX(cols);
        pos(3:4) = [pltW pltH];
        for ros = 1:nRows
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
        end
    end


end
    
    

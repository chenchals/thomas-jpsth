%function [] = corrSpkSatSdfPlot(corrSpkSatSdfFile,pdfOutputDir,savePdfFlag)
%CORRSPKCSATSDFPLOT Summary of this function goes here
%% vars
spkCorrSummaryDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary';
corrSpkSatSdfFile = fullfile(spkCorrSummaryDir,'spkCorrSdfs.mat');
% these are values for column: condition
conditionPairs = {
    {'FastCorrect','AccurateCorrect'};
    {'FastErrorChoice','AccurateErrorChoice'};
    {'FastErrorTiming','AccurateErrorTiming'};
    };
% these are prefixes for sdf/raster column names
sdfColNameSuffix = '_sdfTsMeanStdSem';
epochs = {'Visual','PostSaccade', 'PostReward'};
alignedOnEvents = {'Cue onset','Saccade','Reward'};

%% load data
corrSpkSdfs = load(corrSpkSatSdfFile);
areas = {'SEF','FEF','SC'}; % fields contiong units in corrSpkSatSdfFile
unitInfos = corrSpkSdfs.unitInfo;
filterCriteria = corrSpkSdfs.filterCriteria;
% aggregate
sdfsTbl =table();
for aa = areas
    sdfsTbl = [sdfsTbl;corrSpkSdfs.(aa{1})];     %#ok<*AGROW>
end
clearvars corrSpkSdfs;

%% Plot for each unit in the aggregated table
unitNums = unique(sdfsTbl.unitNum,'stable');
nUnits = numel(unitNums);
for uu = 1:nUnits
    currUnitNum = unitNums(uu);
    currUnitSdfs = sdfsTbl(sdfsTbl.unitNum==currUnitNum,:);
    maxFr = getMaxFr(currUnitSdfs,sdfColNameSuffix);
    currUnitInfo = unitInfos(unitInfos.unitNum==currUnitNum,:);
    % call figure template: plots are column wise, ie, 1st, 2nd,.. cols
    [H_plots,H_Figure] = getPlotHandles();
    % gather sdfs for condition : *Correct, *ErrroChoice, *ErrorTiming
    plotNo = 0;
    for col = 1:numel(epochs)
        epochName = epochs{col};
        alignedOn = alignedOnEvents{col};
        sdfCol = [epochName sdfColNameSuffix];
        
        for ro = 1:numel(conditionPairs)
            currCondPair = conditionPairs{ro};
            condSdfs = currUnitSdfs(ismember(currUnitSdfs.condition,currCondPair),{'condition',sdfCol});
            plotNo = plotNo + 1;
            H_ax = H_plots(plotNo);
            plotSatSdf(H_ax,condSdfs,sdfCol,maxFr);
            xlabel(['Time from ' alignedOn ' (ms)'],'FontWeight','bold');
            % if row = 1 add epoc name at the top
            if ro == 1
                title(epochName);
            end
            % if col == 1, add condition: Correct, ErrorChoice, ErrorTiming
            if col == 1
                annotateCondition(currCondPair);
                ylabel('Spikes/second','FontWeight','bold')
            else
                set(gca,'YColor',[1 1 1]);
            end
        end % for row
    end % for col
    % last column: Paired unit nums with siginif. spike corr for currUnit
    for ro = 1:numel(conditionPairs)
        currCondPair = conditionPairs{ro};
        currUnitSdfs.fastOrAcc = regexprep(currUnitSdfs.condition,'(Correct)|(Error.*)','');
        fastIdx = ismember(currUnitSdfs.condition,currCondPair) & strcmp(currUnitSdfs.fastOrAcc,'Fast');
        accIdx = ismember(currUnitSdfs.condition,currCondPair) & strcmp(currUnitSdfs.fastOrAcc,'Accurate');
        % accurate paired units r_sc <= 0.05
        pairedUnitsTxt = {};
        units = currUnitSdfs.pairedSefUnitNums{accIdx};
        pairedUnitsTxt = [pairedUnitsTxt; getUnitsTxt(units,'SEF','red')];
        units = currUnitSdfs.pairedFefUnitNums{accIdx};
        pairedUnitsTxt = [pairedUnitsTxt; getUnitsTxt(units,'FEF','red')];
        units = currUnitSdfs.pairedScUnitNums{accIdx};
        pairedUnitsTxt = [pairedUnitsTxt; getUnitsTxt(units,'SC','red')];
        
        pairedUnitsTxt =[pairedUnitsTxt; {'  ';'  '}];
        % fast paired units r_sc <= 0.05
        units = currUnitSdfs.pairedSefUnitNums{fastIdx};
        pairedUnitsTxt = [pairedUnitsTxt; getUnitsTxt(units,'SEF','green')];
        units = currUnitSdfs.pairedFefUnitNums{fastIdx};
        pairedUnitsTxt = [pairedUnitsTxt; getUnitsTxt(units,'FEF','green')];
        units = currUnitSdfs.pairedScUnitNums{fastIdx};
        pairedUnitsTxt = [pairedUnitsTxt; getUnitsTxt(units,'SC','green')];
        
        % Position the text
        plotNo = plotNo + 1;
        H_ax = H_plots(plotNo);
        axes(H_ax);
        title({'Paired units';'r_{sc} p<=0.05'});
        set(gca,'XColor',[1 1 1]);
        set(gca,'YColor',[1 1 1]);
        pos = get(gca,'Position');
        annotation('textbox','String',pairedUnitsTxt,'Position',pos,'EdgeColor','none')

    end
    
    %%
    
    delete(gcf)

end % for unit

%%

function [unitsStr] = getUnitsTxt(units,heading,txtColor)
    unitsStr{1,1} = ['\fontsize{10} \color{' txtColor '} \bf' heading ':'];
    if ~isempty(units)
        s = sprintf('%d, ',units);
        unitsStr{2,1} = ['\it   ' s(1:end-2)];
    else
        unitsStr{2,1} = '\it     none ';
    end
end

function [] =  annotateCondition(currCondPair)
    condTxt = regexprep(currCondPair{1},'(Fast)|(Accurate)','');
    xTicks = get(gca,'XTick');
    yTicks = get(gca,'YTick');
    text(min(xTicks)-unique(diff(xTicks)),mean(yTicks),condTxt,...
        'HorizontalAlignment','center','Rotation',90,...
        'FontWeight','bold','FontSize',12,'FontAngle','italic');        
end


function [] = plotSatSdf(H_ax,currSdfTbl,sdfColName,maxFr)

    fx_vecSmooth = @(x,w) smoothdata(x,'movmean',w,'omitnan');
    conds = unique(currSdfTbl.condition,'stable');
    axes(H_ax);
    set(gca,'Box','off');
    pltColor = 'r';
    for c = 1:numel(conds)
        cond = conds{c};
        if contains(cond,'Fast')
            pltColor = 'g';
        end
        currSdf = currSdfTbl.(sdfColName){ismember(currSdfTbl.condition,cond)}(:,1:2);
        plot(currSdf(:,1),fx_vecSmooth(currSdf(:,2),5),'Color',pltColor,'LineWidth',1.5)
        hold on
    end
    ylim([0 maxFr]);
    xTicklabels = get(gca,'XTicklabel');
    xTicklabels{1} = ' ';
    xTicklabels{end} = ' ';
    set(gca,'XTickLabel',xTicklabels);    
    line([0 0],[0 maxFr]);
    set(gca,'Box','off');
end

function [maxFr] = getMaxFr(currUnitSdfs,sdfColNameSuffix)
  colNames = currUnitSdfs.Properties.VariableNames;
  sdfColNames = colNames(~cellfun(@isempty,regexp(colNames,sdfColNameSuffix,'match')));
  allSdfs = cell2mat(cellfun(@(x) cell2mat(currUnitSdfs.(x)),sdfColNames','UniformOutput',false));
  maxFr = ceil(max(allSdfs(:,2)));
end

function [H_plots,H_Figure] = getPlotHandles()
    % 3 condition-rows x 3 epoch-cols-closely-spaced-on-x -> y axis scale only
    % on the leftmost plot
    % plus 1 column for information about pairted units
    % plus 1 row for unit info annotation
     fs = 6;
     if ismac
         fs = 8;
     end
    set(0,'units','pixels');
    set(0,'defaulttextfontsize',fs,...
        'defaulttextfontname','Arial',...
        'defaultaxesfontsize',fs,...
        'defaultaxeslinewidth',0.05);
    margin = 10; %pixels

    ss = [20 20 1500 990];
    FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
    %Main figure window
    H_Figure=figure('Position',FigPos,...
        'color',[1 1 1],'numbertitle','off','renderer','painters',...
        'renderermode','manual','menubar','none',...
        'Tag','H_Figure');
    orient landscape
    set(H_Figure,'units','normalized')

    nRows = 3;
    startOffset = 0.05;
    pltH = 0.22;
    pltWFor3Cols = 0.8;
    gutterx = 0.005;
    guttery = 0.06;
    % conditions SDFs have different lengths of times (x axis).
    % make widths proportional such that time unit ticks are of equal length
    % visual:[-600 400] , postSaccade:[-200 600], postReward:[-100 400]
    pltWProp = [1001 801 501];
    % partition total plot width to the proportion above
    pltWs = pltWProp.*(pltWFor3Cols/sum(pltWProp));

    offsetsX = [0 pltWs(1)+gutterx sum(pltWs(1:2))+gutterx*2] + startOffset;
    offsetsY = 0.90-pltH:-(pltH+guttery):guttery; % for 3 row-starts

    pltCount = 0;
    for col = 1:4
        if col < 4
            pltW = pltWs(col);
            pos(1) = offsetsX(col);
            boxOn = 'on';
       else
            pltW = 1.0 - (startOffset*2 + pltWFor3Cols);
            pos(1) = offsetsX(col-1) + pltWs(col-1) + gutterx*2;
            boxOn = 'off';
        end
        pos(3:4) = [pltW pltH];
        for ro = 1:nRows
            pos(2) = offsetsY(ro);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box',boxOn,...
                'layer','top','Tag',sprintf('H_plot%d',pltCount));
        end
    end
end




% build a bit pattern of connectivity for 
% functional SEF units
% 
% 1. Filter RscTbl for significant connections to SEF units to FEF and SC
%           keep track of positive and negative Rsc
% 2. Get a unique list of significant SEF units and label their VisMovType
% 3. Parse each connection and build a bit pattern for connection
% 4. Bit Pattern for connection: 
% ________________________________________________________________________
%   Bit name   | bit-value  |  Bit meaning
%   FEFconnect | 2^6-1 = 63 | FEF signif. conn.
%   SCconnect  | 2^5-1 = 31 | SC signif. conn.
%   FEFplus1   | 2^4-1 = 15 | FEF Rsc positive, signif. conn. 1 or more
%   FEFminus1  | 2^3-1 = 7  | FEF Rsc negative, signif. conn. 1 or more
%   SCplus1    | 2^2-1 = 3  | SC Rsc positive, signif. conn. 1 or more
%   SCminus1   | 2^1-1 = 1  | SC Rsc negative, signif. conn. 1 or more
% ------------------------------------------------------------------------
rscTbl = getPairRscTableForHeb();

% remove all rows where rhoRaw_150ms is NaN
rscTbl(isnan(rscTbl.rhoRaw_150ms),:) = [];


% code functional types for X
rscTbl.X_visMovType = repmat({'Other'},size(rscTbl,1),1);
rscTbl.visMovTypeOrder = repmat(4,size(rscTbl,1),1);
idx = find(abs(rscTbl.X_visGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Vis'},numel(idx),1);
rscTbl.visMovTypeOrder(idx) = ones(numel(idx),1);
idx = find(abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Mov'},numel(idx),1);
rscTbl.visMovTypeOrder(idx) = 2 + ones(numel(idx),1);
idx = find(abs(rscTbl.X_visGrade) >=2 & abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'VisMov'},numel(idx),1);
rscTbl.visMovTypeOrder(idx) = 1 + ones(numel(idx),1);
% code: rsc sign
rscTbl.rscSign = sign(rscTbl.rhoRaw_150ms);
% code rsc significance : 
% +1 = rsc positive and significant, -1 = rsc negative and significant
rscTbl.signifRsc = zeros(size(rscTbl,1),1);
idx = find(rscTbl.rhoRaw_150ms > 0 & rscTbl.pvalRaw_150ms <= 0.01);
rscTbl.signifRsc(idx) = ones(numel(idx,1));
idx = find(rscTbl.rhoRaw_150ms < 0 & rscTbl.pvalRaw_150ms <= 0.01);
rscTbl.signifRsc(idx) = ones(numel(idx,1));
sigRscTbl = rscTbl(rscTbl.signifRsc == 1,:);

% filter rsc table for for SEF to FEF/SC connections only
rscTbl(~ismember(rscTbl.X_area,{'SEF'}),:) = [];
rscTbl(ismember(rscTbl.Y_area,{'SEF'}),:) = [];
sigRscTbl(~ismember(sigRscTbl.X_area,{'SEF'}),:) = [];
sigRscTbl(ismember(sigRscTbl.Y_area,{'SEF'}),:) = [];



%% get counts and distributions for paper only SEF to FEF/SC connections
r = rscTbl(ismember(rscTbl.epoch,'PostSaccade'),:);
s = sigRscTbl(ismember(sigRscTbl.epoch,'PostSaccade'),:);
crossAreaTblPs = table();
crossAreaTblPs.sefType = {'Vis','VisMov','Mov','Other'}';
crossAreaTblPs.X_visMovType = {'Vis','VisMov','Mov','Other'}';

t = unique(r(:,{'X_visMovType','X_unitNum'}));
temp = grpstats(t,'X_visMovType');
temp = join(crossAreaTblPs,temp,'Keys','X_visMovType');
crossAreaTblPs.sefCount = temp.GroupCount;
t = unique(s(:,{'X_visMovType','X_unitNum'}));
temp = grpstats(t,'X_visMovType');
temp = join(crossAreaTblPs,temp,'Keys','X_visMovType');
crossAreaTblPs.sefCountSig = temp.GroupCount;

%% pair counts for FEF/SC by sef type
for z = {'FEF','SC'}
    a = z{1};
    f = unique(r(ismember(r.Y_area,a),{'X_visMovType','X_unitNum','Y_unitNum'}));
    fs = unique(s(ismember(s.Y_area,a),{'X_visMovType','X_unitNum','Y_unitNum'}));
    temp = grpstats(f,'X_visMovType');
    tempSum = sum(temp.GroupCount);
    vn = sprintf('%sPairs_%d_%d',a,numel(unique(f.Y_unitNum)),tempSum);
    temp.(vn) = temp.GroupCount;
    crossAreaTblPs = join(crossAreaTblPs,temp,'Keys','X_visMovType','RightVariables',vn);
    temp = grpstats(fs,'X_visMovType');
    tempSum = sum(temp.GroupCount);
    vn = sprintf('%sPairsSig_%d_%d',a,numel(unique(f.Y_unitNum)),tempSum);
    temp.(vn) = temp.GroupCount;
    crossAreaTblPs = join(crossAreaTblPs,temp,'Keys','X_visMovType','RightVariables',vn);
end
crossAreaTblPs.sefType{5} = 'All';
crossAreaTblPs.X_visMovType{5} = 'All';
crossAreaTblPs{5,3:end} = sum(crossAreaTblPs{:,3:end},1);
crossAreaTblPs.X_visMovType = [];
writetable(crossAreaTblPs,'crossAreaTable.csv')

%% build the word pattern table for significant connections

connPatternTbl = getConnWordPattern(sigRscTbl);
% uniqDecToBinPattern = unique(wordPattern(:,{'connWord6bitDec', 'connWord4bitDec','connWord6bit','connWord4bit'}));
funcWordPattern = getFuncWordPattern(connPatternTbl);

%% permute functional type labels on the connPatternTbl to simulate random functional type associations
sefUnitType = unique(sigRscTbl(:,{'X_unitNum','X_visMovType'}));
temp = sefUnitType.X_unitNum;
nPerms = 10;
randFuncWordPatternPerms = cell(nPerms,1);
for np = 1:nPerms
    fprintf('premutation #%d\n',np);
    % permute unit labels
    randSefType = sefUnitType.X_visMovType(randperm(numel(temp)));
    fn = ['rand_' num2str(np,'%d')];
    sefUnitType.(fn) = randSefType;
    % assign this random func type to table variable
    rndRscTbl = sigRscTbl;
    for u = 1:numel(temp)
        uNum = temp(u);
        nType = sefUnitType.(fn)(sefUnitType.X_unitNum == uNum);
        idx = rndRscTbl.X_unitNum == uNum;
        rndRscTbl.X_visMovType(idx) = repmat(nType,sum(idx),1);
    end
    Z = getConnWordPattern(rndRscTbl);
    randFuncWordPatternPerms{np} = getFuncWordPattern(Z);
end
% Create mean_GroupCount for the randFuncWordPattern 
randFuncWordPatternAll = table();
for ii = 1:size(randFuncWordPatternPerms,1)
    randFuncWordPatternAll = [randFuncWordPatternAll;randFuncWordPatternPerms{1}];
end
useVars = {'outcome','epoch','satCondition','sefVisMovType','connWord6bitDec','connWord4bitDec','connWord6bit','connWord4bit','GroupCount'};

randFuncWordPattern = grpstats(randFuncWordPatternAll(:,useVars),useVars(1:end-1),'mean');
randFuncWordPattern.Properties.RowNames = {};
colNames = randFuncWordPattern.Properties.VariableNames;
colNames = regexprep(colNames,{'^GroupCount$','^mean_GroupCount$'},{'nPerms','GroupCount'});
randFuncWordPattern.Properties.VariableNames = colNames;

%% now we have 
% funcWordPattern and randFuncWordPattern also sefUnitType for what random
% combinations were used for randFuncWordPattern
% remove all vars except these
clearvars -except funcWordPattern randFuncWordPattern sefUnitType
% create all possible patterns There will be 15 unique patters
% 
uniqPatternTbl = table();
for ii = 1:15
    pat = table();
    p = dec2bin(ii,4);
    pat.patIdx = ii;
    pat.connWord4bitDec = ii;
    pat.connWord4bit = {p};
    pat.FEFplus1 = p(1);
    pat.FEFminus1 = p(2);
    pat.SCplus1 = p(3);
    pat.SCminus1 = p(4);
    uniqPatternTbl =[uniqPatternTbl;pat];
end


patColNames = strcat('pat_',num2str(uniqPatternTbl.connWord4bitDec,'%02d_'),uniqPatternTbl.connWord4bit)';
featureColNames = {'outcome','epoch','satCondition','sefVisMovType'};
permuteColNames = {'dataSrc','nPrems'};
useCols = [permuteColNames featureColNames patColNames];
patternSummTbl = table();

tbls = {funcWordPattern,randFuncWordPattern};
for ii = 1:numel(tbls)
    dataSrc = 'actual';
    nPerms = 0;
    dat = tbls{ii};
    if ii == 2
        dataSrc = 'randPerm';
        nPerms = unique(dat.nPerms);
    end
    % build table row for each entry in the dat table
    nRows = size(dat,1);
    for nr = 1:nRows
        datRo = dat(nr,:);
        outcome = datRo.outcome;
        epoch = datRo.epoch;
        satCondition = datRo.satCondition;
        sefVisMovType = datRo.sefVisMovType;
        groupCount = datRo.GroupCount;
        pattern4colName = char(strcat('pat_', num2str(datRo.connWord4bitDec,'%02d_'), datRo.connWord4bit));
        
        % find idx into table to aggregate across patterns
        idx = 0;
        if ~isempty(patternSummTbl)
        idx = find(ismember(patternSummTbl.dataSrc,dataSrc) & patternSummTbl.nPerms == nPerms...
            & ismember(patternSummTbl.outcome,outcome) & ismember(patternSummTbl.epoch,epoch)...
            & ismember(patternSummTbl.satCondition,satCondition)...
            & ismember(patternSummTbl.sefVisMovType,sefVisMovType));
        end
        if idx > 0
            patternSummTbl.(pattern4colName)(idx) = groupCount;
        else
            p = table();
            p.dataSrc{1} = dataSrc;
            p.nPerms = nPerms;
            p.outcome = outcome;
            p.epoch = epoch;
            p.satCondition = satCondition;
            p.sefVisMovType = sefVisMovType;
            p = [p array2table(zeros(1,numel(patColNames)),'VariableNames',patColNames)];
            p.(pattern4colName) = groupCount;
            patternSummTbl = [patternSummTbl; p];
        end
    end
end

% sum total count across rows
patternSummTbl.nAllPatterns = sum(patternSummTbl{:,patColNames},2);

%clearvars -except funcWordPattern randFuncWordPattern sefUnitType uniqPatternTbl patternSummTbl

fastCorrectPostScacc = patternSummTbl.dataSrc



%% show a heatmap of patterns
% {'actual',0,'ErrorTiming','PostSaccade','Accurate'}
dataSrcs = {'actual','randPerm'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
epochs = {'PostSaccade'};
satConds = {'Accurate','Fast'};
allDat = patternSummTbl;
colNames = allDat.Properties.VariableNames;
funcLabelsAll = {'Vis','VisMov','Mov','Other'};
datColnames = colNames(~cellfun(@isempty, regexp(colNames,'pat_.*','match')));
cLims = [1 max(max(allDat{:,datColnames}))];
for ds = 1:numel(dataSrcs)
    h_fig = newFigure();
    pltNo = 0;
    dataSrc = dataSrcs{ds};
    annotation('textbox','Position',[0.05 0.95,0.9,0.05],'String',dataSrc,'FontSize',16,'FitBoxToText','off','EdgeColor','none','FontWeight','bold')
    idxDs = ismember(allDat.dataSrc,dataSrc);
    rowSum = {};
    colSum = {};
    for ep = 1:numel(epochs)
        epoch = epochs{ep};
        idxEpoch = ismember(allDat.epoch,epoch);
        for oc = 1:numel(outcomes)
            outcome = outcomes{oc};
            outcomeEpoch =  [outcome '-' epoch];
            %annotation('textbox','Position',[0.15 0.95,0.8,0.05],'String',outcomeEpoch,'FontSize',14,'FitBoxToText','off','EdgeColor','none','FontWeight','bold')
            idxOutcome = ismember(allDat.outcome,outcome);
            dat = allDat(idxDs & idxEpoch & idxOutcome,:);
            idxFast = ismember(dat.satCondition,'Fast');
            idxAcc = ismember(dat.satCondition,'Accurate');
            % for Fast
            [heatmapDat,funcLabels,patLabels] = ...
                getHeatmapDat(dat(idxFast,:),funcLabelsAll,datColnames);
            pltNo = pltNo + 1;
            h_axis = subplot(1,6,pltNo);
            fx_plotHeatmap(h_axis,heatmapDat,funcLabels,patLabels,cLims) 
            title(['Fast-' outcomeEpoch])
            % for Accurate
            [heatmapDat,funcLabels,patLabels] = ...
                getHeatmapDat(dat(idxAcc,:),funcLabelsAll,datColnames);
            pltNo = pltNo + 1;
            h_axis = subplot(1,6,pltNo);
            fx_plotHeatmap(h_axis,heatmapDat,funcLabels,patLabels,cLims);
            title(['Accurate-' outcomeEpoch])
                       
        end
    end
    %saveFigPdf(['connectionPattern-' dataSrc '.pdf'])
end

%% Create graph and subgraphs from the patternSummTbl
% Extract fromUnitGrp and toUnitGrp tags:
% 1. code fromUnitGrp
temp = patternSummTbl;
nR = size(temp,1);
preFixCols = {'dataSrc','outcome','epoch','satCondition','sefVisMovType'};
grpPreFix = arrayfun(@(x) sprintf('%s-%s-%s%s',char(join(temp{x,preFixCols},'-')),...
                          regexprep(temp{1,'sefVisMovType'},'[a-z]',''),num2str(x,'-%02d')),...
                          (1:nR)','UniformOutput',false);

% 2. 
      

% 3. 

%%
function [heatmapDat,funcLabels,patLabels] = getHeatmapDat(inDatTbl,funcLabelsAll,datColnames)
    missingFuncs = setdiff(funcLabelsAll,inDatTbl.sefVisMovType);
    for mf = 1:numel(missingFuncs)
        inDatTbl.sefVisMovType{size(inDatTbl,1)+1} = missingFuncs{mf};
    end
    % sort funcs
    for ft = 1:numel(funcLabelsAll)
        idx = ismember(inDatTbl.sefVisMovType,funcLabelsAll{ft});
        inDatTbl.sortOrder(idx) = ft;
    end
    inDatTbl = sortrows(inDatTbl,'sortOrder');
    heatmapDat = inDatTbl{:,datColnames};
    funcLabels = inDatTbl.sefVisMovType;
    patLabels = regexprep(datColnames,'pat_\d+_','');
    
end

function [h_heatmap,rowSum,colSum] = fx_plotHeatmap(h_axis,dat,fxLabels,patLabels,cLims)  
    cData = dat;
    cData(cData == 0) = NaN;
    rowSum = nansum(cData,2);
    colSum = nansum(cData,1)';
    pLabels = strcat(patLabels',num2str(colSum,'(%d)'));   
    funcLabels = strcat(fxLabels,num2str(rowSum,' (%d)'));    
    
    axes(h_axis)
    h_heatmap = heatmap(funcLabels,pLabels,cData','ColorLimits',cLims);

    colormap('cool')
    h_heatmap.MissingDataColor = [0.95 0.95 0.95];
    h_heatmap.MissingDataLabel = 'No Conn.';
    h_heatmap.ColorbarVisible = 'off';
end



%% Get word pattern tables for original input Rsc table
function [connWordPattern] = getConnWordPattern(inRscTbl)
    % filter table for outcomeEpoch
    epochs = {'PostSaccade'};
    outcomes = {'Correct','ErrorChoice','ErrorTiming'};
    satConds = {'Fast','Accurate'};

    connWordPattern = table();
    for ep = 1:numel(epochs)
        epoch = epochs{ep};
        for oc = 1:numel(outcomes)
            outcome = outcomes{oc};
            rscs = inRscTbl(ismember(inRscTbl.epoch,epoch) & ismember(inRscTbl.outcome,outcome),:);
            for sc = 1:numel(satConds)
                satCond = satConds{sc};
                dat = rscs(ismember(rscs.satCondition,satCond),:);
                sefUnitNums = unique(dat.X_unitNum);
                for u = 1:numel(sefUnitNums)
                    unitNum = sefUnitNums(u);
                    unitPairs = dat(dat.X_unitNum == unitNum,:);
                    plusIdx = unitPairs.rscSign == 1;
                    minusIdx = unitPairs.rscSign == -1;
                    % Count FEF plus/minus pairs
                    fefNconnPlus = sum(ismember(unitPairs.Y_area,'FEF') & plusIdx);
                    fefNconnMinus = sum(ismember(unitPairs.Y_area,'FEF') & minusIdx);
                    % Count SC plus/minus pairs
                    scNconnPlus = sum(ismember(unitPairs.Y_area,'SC') & plusIdx);
                    scNconnMinus = sum(ismember(unitPairs.Y_area,'SC') & minusIdx);
                    % build a table row of output
                    wTbl = table();
                    wTbl.outcome{1} = outcome;
                    wTbl.epoch{1} = epoch;
                    wTbl.satCondition{1} = satCond;
                    wTbl.sefVisMovType = unique(unitPairs.X_visMovType);
                    wTbl.sefUnitNum = unitNum;
                    wTbl.nSignifPairs = size(unitPairs,1);
                    wTbl.nSignifPlusRsc = sum(plusIdx);
                    wTbl.nSignifMinusRsc = sum(minusIdx);
                    % parse the pairs for this unit into a bit pattern word
                    % initialize bit pattern
                    wTbl.FEFconnected = 0;
                    wTbl.SCconnected = 0;
                    wTbl.FEFplus1 = 0;
                    wTbl.FEFminus1 = 0;
                    wTbl.SCplus1 = 0;
                    wTbl.SCminus1 = 0;
                    % for FEF pairs
                    if fefNconnPlus > 0
                        wTbl.FEFconnected = 1;
                        wTbl.FEFplus1 = 1;
                    end
                    if fefNconnMinus > 0
                        wTbl.FEFconnected = 1;
                        wTbl.FEFminus1 = 1;
                    end
                    % for SC pairs
                    if scNconnPlus > 0
                        wTbl.SCconnected = 1;
                        wTbl.SCplus1 = 1;
                    end
                    if scNconnMinus > 0
                        wTbl.SCconnected = 1;
                        wTbl.SCminus1 = 1;
                    end
                    connWordPattern = [connWordPattern; wTbl]; %#ok<*AGROW>
                    % done for this unit
                end
            end
        end
    end

    % add column for 6 bit word
    nRows = size(connWordPattern,1);
    colsFor6bitWord = {'FEFconnected','SCconnected','FEFplus1','FEFminus1','SCplus1','SCminus1'};
    colsFor4bitWord = {'FEFplus1','FEFminus1','SCplus1','SCminus1'};
    connWordPattern.connWord6bit = arrayfun(@(x) sprintf('%d',connWordPattern{x,colsFor6bitWord}),(1:nRows)','UniformOutput',false);
    connWordPattern.connWord4bit = arrayfun(@(x) sprintf('%d',connWordPattern{x,colsFor4bitWord}),(1:nRows)','UniformOutput',false);
    % convert bit pattern to decimal
    connWordPattern.connWord6bitDec = arrayfun(@(x) bin2dec(x{1}),connWordPattern.connWord6bit);
    connWordPattern.connWord4bitDec = arrayfun(@(x) bin2dec(x{1}),connWordPattern.connWord4bit);

end



%% Count different patterns by SEF functional types
function [funcWordTbl] = getFuncWordPattern(connPatternTbl)
    epochs = {'PostSaccade'};
    outcomes = {'Correct','ErrorChoice','ErrorTiming'};
    satConds = {'Fast','Accurate'};
    useVars = {'outcome','epoch','satCondition','sefVisMovType','connWord6bitDec','connWord4bitDec','connWord6bit','connWord4bit'};
    grpVars = useVars; %same order
    funcWordTbl = table();
    for oc = 1:numel(outcomes)
        outcome = outcomes{oc};
        idxOutcome = ismember(connPatternTbl.outcome,outcome);
        for ep = 1:numel(epochs)
            epoch = epochs{ep};
            idxEpoch = ismember(connPatternTbl.epoch,epoch);
            for sc = 1:numel(satConds)
                satCond = satConds{sc};
                idx =  ismember(connPatternTbl.satCondition,satCond) & idxOutcome & idxEpoch;
                temp = connPatternTbl(idx,:);
                res = grpstats(temp(:,useVars),grpVars,{'counts'});
                res.Properties.RowNames = {};
                funcWordTbl = [funcWordTbl;res];
            end
        end
    end
end



              


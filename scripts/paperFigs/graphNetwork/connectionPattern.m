% Do the following:
% 1. Aggregate units as follows to build a matrix
%      SEF: VIS, VISMOV, MOV, OTHER categories
%      FEF: FEF_1, FEF_2
%      SC: SC_1, SC_2
%      
% 2. Increment significant connection counts

rscTbl = getPairRscTableForHeb();

% remove all rows where rhoRaw_150ms is NaN
rscTbl(isnan(rscTbl.rhoRaw_150ms),:) = [];


% filter rsc table for for SEF to FEF/SC connections only
rscTbl(~ismember(rscTbl.X_area,{'SEF'}),:) = [];
rscTbl(ismember(rscTbl.Y_area,{'SEF'}),:) = [];

% code functional types for X
rscTbl.X_visMovType = repmat({'Other'},size(rscTbl,1),1);
idx = find(abs(rscTbl.X_visGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Vis'},numel(idx),1);
idx = find(abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Mov'},numel(idx),1);
idx = find(abs(rscTbl.X_visGrade) >=2 & abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'VisMov'},numel(idx),1);
% code: rsc sign
rscTbl.rscSign = sign(rscTbl.rhoRaw_150ms);
% code rsc significance : 
% +1 = rsc positive and significant, -1 = rsc negative and significant
rscTbl.signifRsc = zeros(size(rscTbl,1),1);
idx = find(rscTbl.rhoRaw_150ms > 0 & rscTbl.pvalRaw_150ms <= 0.01);
rscTbl.signifRsc(idx) = ones(numel(idx,1));
idx = find(rscTbl.rhoRaw_150ms < 0 & rscTbl.pvalRaw_150ms <= 0.01);
rscTbl.signifRsc(idx) = -1*ones(numel(idx,1));

%% count number of significant connections to build a heatmap matrix
% 1. colNames =
%           {'satCondition','rscSign','sefConnType','SEFvisMoveType','nPairs','nSignif',...
%            'FEF_Grp1','FEF_Grp2','SC_Grp1','SC_Grp2',... % split nConnections as 1 + (n-1)
%            'FEF_1','FEF_2','SC_1','SC_2'... % split nConnections n/2
%            };
%
%           satCondition = {'Fast','Accurate'}
%           rscSign = {'all','positive','negative'}
%           sefConnType = {'single','multiple'} does a single sef unit
%                         have one or more signif connections?
%           SEFvisMoveType = {'Vis','VisMov','Mov','Other'}
%           nPairs = number of pairs (signif and non-signif)
%           nSignif = number of pairs that are significant
% 2. nows = 2 (satCondition) 
%         * 4 (SEFvisMoveType) 
%         * 3 (rscSign)
%         * 2 (sefConnType) 
%       = 48 rows for each outcome epoch combination
% 3. updating matrix-entry for a row: 
%   Case1: Partial Connection count:  
%      Example: current SEF unit is Vis then
%      a. number of significant connections nc == 1
%         i. if connection is to FEF then 
%                update FEF_1 = FEF_1 + 0.5, FEF_2 = FEF_2 + 0.5
%         ii. if connection is to SC then 
%                update SC_1 = SC_1 + 0.5, SC_2 = SC_2 + 0.5
%      b. number of significant connections nc > 1
%         i. Split nc connections by area: nc_fef and nc_sc
%         ii. update columns FEF_1 = FEF_1 + nc_fef/2, FEF_2 = FEF_2 + nc_fef/2
%         ii. update columns SC_1 = SC_1 + nc_sc/2, SC_2 = SC_2 + nc_sc/2
%   
%   Case2: Whole connection count:
%      Example: current SEF unit is Vis then
%      a. number of significant connections nc == 1
%         i. if connection is to FEF then 
%                update FEF_grp1 = FEF_grp1 + nc
%         ii. if connection is to SC then 
%                update SC_grp1 = SC_grp1 + nc
%      b. number of significant connections nc > 1
%         i. Split nc connections by area: nc_fef and nc_sc
%         ii. update columns if nc_fef > 0 then 
%                             FEF_Grp1 = FEF_Grp1 + 1,
%                             FEF_Grp2 = if (nc_fef - 1) > 0, then FEF_Grp2 +  (nc_fef - 1)
%         ii. update columns if nc_sc > 0 then
%                             SC_Grp1 = SC_Grp1 + 1,
%                             SC_Grp2 = if (nc_sc - 1) > 0, then SC_Grp2 +  (nc_sc - 1)
%
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
epochs = {'PostSaccade'};
satConds = {'Fast','Accurate'};
rscSigns = {'all','positive','negative'};
sefConnTypes = {'single','multiple'};
SEFVisMovType = {'Vis','VisMov','Mov','Other'};
colValsStr = {'satCondition','rscSign','sefConnType','SEFvisMoveType'};
colValsNum = {'nPairs','nSignif','FEF_Grp1','FEF_Grp2','SC_Grp1','SC_Grp2','FEF_1','FEF_2','SC_1','SC_2'};
% Create an empty connection table 
% for initiatizing tables in the loop below
emptyConnTbl = table();
for rscSign = rscSigns
    for connType = sefConnTypes
        for sCond = satConds
            temp = table();
            nr = numel(SEFVisMovType);
            temp.satCondition = repmat(sCond,nr,1);
            temp.rscSign = repmat(rscSign,nr,1);
            temp.sefConnType = repmat(connType,nr,1);
            temp.SEFvisMoveType = SEFVisMovType';
            for cn = colValsNum
                temp.(cn{1}) = zeros(nr,1);
            end
            emptyConnTbl = [emptyConnTbl;temp]; %#ok<*AGROW>
        end
    end
end
%%
outConnMat = struct();

for oc = 1:numel(outcomes)
   outcome = outcomes{oc};
    rscs = rscTbl(ismember(rscTbl.outcome,outcome),:);
   for ep = 1:numel(epochs)
       epoch = epochs{ep};
       rscs = rscs(ismember(rscs.epoch,epoch),:);
       outFieldname = [outcome epoch];
       outConn = emptyConnTbl;
       for sc = 1:numel(satConds)
           satCond = satConds{sc};
           fastAccuPairs = rscs(ismember(rscs.satCondition,satCond),:);
           % Process connections for each unit
           sefUnits = unique(rscs.X_unitNum);
           for u = 1:numel(sefUnits)
               unitNum = sefUnits(u);
               currUnitPairs = fastAccuPairs(fastAccuPairs.X_unitNum == unitNum,:);
               visMovType = char(unique(currUnitPairs.X_visMovType));
                   % process for all Rsc do not take sign into account
                   for rsig = 1:numel(rscSigns)
                       rscSign = rscSigns{rsig};                       
                       switch rscSign
                           case 'positive'
                               tempPairs = currUnitPairs(currUnitPairs.rscSign == 1,:);                              
                           case 'negative'
                               tempPairs = currUnitPairs(currUnitPairs.rscSign == -1,:);                                                       
                           case 'all'
                               tempPairs = currUnitPairs(abs(currUnitPairs.rscSign) == 1,:);
                       end
                       
                       if isempty(tempPairs)
                           continue;
                       end
                       % assume there is only 1 connection for this unit
                       nPairs = size(tempPairs,1);
                       % retrieve all row index of output table {single'
                       % and 'multiple'
                       idx = find(ismember(outConn.satCondition,satCond) ...
                                  & ismember(outConn.SEFvisMoveType,visMovType) ...
                                  & ismember(outConn.rscSign,rscSign));
                       % update values
                       outConn.nPairs(idx) = outConn.nPairs(idx) + nPairs;
                         
                       % Prune to significant connections
                       tempPairs(tempPairs.signifRsc == 0,:) = [];                       
                       if isempty(tempPairs)
                           continue; % no significant pairs
                       end
                       nSignif = size(tempPairs,1);
                       sefConnType = 'single';
                       if nSignif > 1
                           sefConnType = 'multiple';
                       end
                       % update output table
                       % retrieve row index of output table
                       idx = find(ismember(outConn.satCondition,satCond) ...
                                  & ismember(outConn.SEFvisMoveType,visMovType) ...
                                  & ismember(outConn.rscSign,rscSign)...
                                  & ismember(outConn.sefConnType,sefConnType));
                       
                       % update nSignif values in the retrieved row
                       outConn.nSignif(idx) = outConn.nSignif(idx) + nSignif;
                       % update FEF_1 and FEF_2 values
                       nc = numel(find(ismember(tempPairs.Y_area,'FEF')));
                       % case1:
                       outConn.FEF_1(idx) = outConn.FEF_1(idx) + nc/2;
                       outConn.FEF_2(idx) = outConn.FEF_2(idx) + nc/2;
                       if nc > 0
                           % case2:
                           outConn.FEF_Grp1(idx) = outConn.FEF_Grp1(idx) + 1;
                           nc = nc - 1;
                           nc = (abs(nc) + nc)/2; % nc is 0 or positive
                           outConn.FEF_Grp2(idx) = outConn.FEF_Grp2(idx) + nc;
                       end
                       % update SC_1 and SC_2 values
                       % case1:
                       nc = numel(find(ismember(tempPairs.Y_area,'SC')));
                       outConn.SC_1(idx) = outConn.SC_1(idx) + nc/2;
                       outConn.SC_2(idx) = outConn.SC_2(idx) + nc/2;
                       % case2:
                       if nc > 0
                           outConn.SC_Grp1(idx) = outConn.SC_Grp1(idx) + 1;
                           nc = nc - 1;
                           nc = (abs(nc) + nc)/2; % nc is 0 or positive
                           outConn.SC_Grp2(idx) = outConn.SC_Grp2(idx) + nc;
                       end
                       
                  end
 
               % next unit num
           end
           outConnMat.(outFieldname) = outConn;
 
       end
   end
end


%% Create heatmaps? for single and multiple connections
% redifine stuff for creating heat maps
fns = fieldnames(outConnMat); % CorrectPostSaccade, ErrorChoicePostSaccade, ErrorTimingPostSaccade
sefConnTypes = {'single','multiple'};
satCons = {'Fast','Accurate'};
rscSigns = {'all','positive','negative'};
heatmapXCols = {'FEF_Grp1','FEF_Grp2','SC_Grp1','SC_Grp2'};
heatmapYCols = {'SEFvisMoveType'};
filterCols = {'satCondition','rscSign','sefConnType'};
plotCols = [filterCols heatmapYCols heatmapXCols];
fx_dat_filter = @(t,satCond,sefConnType,rscSign) t(ismember(t.sefConnType,sefConnType) ...
                                     & ismember(t.rscSign,rscSign)...
                                     & ismember(t.satCondition,satCond)...
                                    ,:);
                                
% find min/max values for color...
tempVals = [];
for oe =1 :numel(fns)
    tempVals = [tempVals;outConnMat.(fns{oe}){:,heatmapXCols}];
end
tempVals = tempVals(:);
cLims = [min(tempVals),max(tempVals)];

for oe = 1:numel(fns)
    outcomeEpoch = fns{oe};
    tempTbl = outConnMat.(outcomeEpoch)(:,plotCols);
    % plot single connection FAST/ACCU for all, positive, negative Rsc
    f = newFigure();
    annotation('textbox','Position',[0.05 0.95,0.9,0.05],'String',outcomeEpoch,'FontSize',14,'FitBoxToText','off','EdgeColor','none')
    pltNo = 0;
    for ro = 1:3
        rscSign = rscSigns{ro};
        for colGrp = 1:2
            sefConnType = sefConnTypes{colGrp};
            for fa = 1:2
                satCond = satConds{fa};
                dat = fx_dat_filter(tempTbl,satCond,sefConnType,rscSign);
                pltNo = pltNo + 1;
                h_axes(pltNo) = subplot(3,4,pltNo);
                h_heatmaps(pltNo) = fx_plotHeatmap(h_axes(pltNo),dat,heatmapXCols,cLims);
                title([satCond '-' rscSign '-' sefConnType])
            end
        end
    end  
    saveFigPdf([outcomeEpoch '.pdf']);
    delete(f)
end


%%
function [h_heatmap] = fx_plotHeatmap(h_axis,dat,heatmapXCols,cLims) 
    cData = dat{:,heatmapXCols};
    xLabel = regexprep(heatmapXCols','_Grp','');
    yLabel = dat.SEFvisMoveType;
    axes(h_axis)
    h_heatmap = heatmap(xLabel,yLabel,cData,'ColorLimits',cLims);
    colormap('cool')
end


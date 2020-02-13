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
rscTbl.signifRsc(idx) = ones(numel(idx,1));

%% count number of significant connections to build a heatmap matrix
% 1. colNames =
%           {'satCondition','rscSign','sefConnType','SEFvisMoveType','nPairs','nSignif',...
%            'FEF_Grp1','FEF_Grp2','SC_Grp1','SC_Grp2',... % split nConnections as 1 + (n-1)
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
satConds = {'Accurate','Fast'};
rscSigns = {'positive','negative','all'};
sefConnTypes = {'multiple','single'};
sefVisMovTypes = {'Vis','VisMov','Mov','Other'}';% row values
colValsStr = {'satCondition','rscSign','sefConnType','SEFvisMoveType'};
colValsNum = {'nPairs','nSignif','FEF_Grp1','FEF_Grp2','SC_Grp1','SC_Grp2'};
% Create an empty connection table 
% for initiatizing tables in the loop below
emptyConnTbl = table();
for rscSign = rscSigns
    for connType = sefConnTypes
        for sCond = satConds
            temp = table();
            nr = numel(sefVisMovTypes);
            temp.satCondition = repmat(sCond,nr,1);
            temp.rscSign = repmat(rscSign,nr,1);
            temp.sefConnType = repmat(connType,nr,1);
            temp.sefVisMovType = sefVisMovTypes;
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
               unitPairs = fastAccuPairs(fastAccuPairs.X_unitNum == unitNum,:);
               % since each unit has a mutually exclusive functional type
               % no for loop for processing
               visMovType = char(unique(unitPairs.X_visMovType));
               % rsc signs +, -, all; use for loop...
               for rsig = 1:numel(rscSigns)
                   rscSign = rscSigns{rsig};
                   switch rscSign
                       case 'positive'
                           signPairs = unitPairs(unitPairs.rscSign == 1,:);
                       case 'negative'
                           signPairs = unitPairs(unitPairs.rscSign == -1,:);
                       case 'all'
                           signPairs = unitPairs;
                   end
                   
                   if ~isempty(signPairs)
                       % all significant and non-significant pairs
                       nPairs = size(signPairs,1);
                       % output table row idx for  {'single','multiple'}
                       outTblRoIdx = find(strcmp(outConn.satCondition,satCond) ...
                           & strcmp(outConn.sefVisMovType,visMovType) ...
                           & strcmp(outConn.rscSign,rscSign));
                       % update values for signed pairs
                       outConn.nPairs(outTblRoIdx) = outConn.nPairs(outTblRoIdx) + nPairs;
                       
                       % process only significant pairs
                       signifPairs = signPairs(signPairs.signifRsc == 1,:);
                       nSignif = size(signifPairs,1);
 
                       if nSignif > 0
                           if nSignif > 1
                               sefConnType = 'multiple';
                           else
                               sefConnType = 'single';
                           end
                           % update output table
                           % retrieve row index of output table
                           outTblRoIdx = find(strcmp(outConn.satCondition,satCond) ...
                               & strcmp(outConn.sefVisMovType,visMovType) ...
                               & strcmp(outConn.rscSign,rscSign)...
                               & strcmp(outConn.sefConnType,sefConnType));
                           
                           % update nSignif values in the retrieved row
                           outConn.nSignif(outTblRoIdx) = outConn.nSignif(outTblRoIdx) + nSignif;
                           % update FEF_1 and FEF_2 values
                           nc = sum(strcmp(signifPairs.Y_area,'FEF'));
                           if nc > 0
                               % case2:
                               outConn.FEF_Grp1(outTblRoIdx) = outConn.FEF_Grp1(outTblRoIdx) + 1;
                               if (nc > 1)
                                   outConn.FEF_Grp2(outTblRoIdx) = outConn.FEF_Grp2(outTblRoIdx) + nc - 1;
                               end
                           end
                           % update SC_1 and SC_2 values
                           nc = sum(strcmp(signifPairs.Y_area,'SC'));
                           if nc > 0
                               % case2:
                               outConn.SC_Grp1(outTblRoIdx) = outConn.SC_Grp1(outTblRoIdx) + 1;
                               if (nc > 1)
                                   outConn.SC_Grp2(outTblRoIdx) = outConn.SC_Grp2(outTblRoIdx) + nc - 1;
                               end
                           end
                           %if ~isempty(signifPairs)
                       end
                       % ~isempty(signPairs)
                   end
                   % next rscSign
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
satConds = {'Fast','Accurate'};
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
cLims = [1,max(tempVals)];

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
    %saveFigPdf([outcomeEpoch '.pdf']);
    %delete(f)
end


%%
function [h_heatmap] = fx_plotHeatmap(h_axis,dat,heatmapXCols,cLims) 
    cData = dat{:,heatmapXCols};
    cData(cData == 0) = NaN;
    xLabel = regexprep(heatmapXCols','_Grp','');
    yLabel = dat.sefVisMovType;
    axes(h_axis)
    h_heatmap = heatmap(xLabel,yLabel,cData,'ColorLimits',cLims);
    colormap('cool')
    h_heatmap.MissingDataColor = [0.95 0.95 0.95];
    h_heatmap.MissingDataLabel = 'No Conn.';
    h_heatmap.ColorbarVisible = 'off';
end


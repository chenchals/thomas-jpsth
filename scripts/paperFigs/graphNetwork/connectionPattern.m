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
%           {'satCondition','rscSign','sefConnType','SEFvisMoveType','nPairs','nSignif','FEF_1','FEF_2','SC_1','SC_2'}
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
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
epochs = {'PostSaccade'};
satConds = {'Fast','Accurate'};
rscSigns = {'all','positive','negative'};
sefConnTypes = {'single','multiple'};
SEFVisMovType = {'Vis','VisMov','Mov','Other'};
colNames = {'satCondition','rscSign','sefConnType','SEFvisMoveType','nPairs','nSignif','FEF_1','FEF_2','SC_1','SC_2'};
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
            for cn = {'nPairs','nSignif','FEF_1','FEF_2','SC_1','SC_2'}
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
           rscs = rscs(ismember(rscs.satCondition,satCond),:);
           % Process connections for each unit
           sefUnits = unique(rscs.X_unitNum);
           for u = 1:numel(sefUnits)
               unitNum = sefUnits(u);
               currPairs = rscs(rscs.X_unitNum == unitNum,:);
               visMovType = char(unique(currPairs.X_visMovType));
                   % process for all Rsc do not take sign into account
                   for rsig = 1:numel(rscSigns)
                       rscSign = rscSigns{rsig};                       
                       switch rscSign
                           case 'positive'
                               tempPairs = currPairs(currPairs.rscSign == 1,:);                              
                           case 'negative'
                               tempPairs = currPairs(currPairs.rscSign == -1,:);                                                       
                           case 'all'
                               tempPairs = currPairs(abs(currPairs.rscSign) == 1,:);
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
                       outConn.FEF_1(idx) = outConn.FEF_1(idx) + nc/2;
                       outConn.FEF_2(idx) = outConn.FEF_2(idx) + nc/2;
                       % update SC_1 and SC_2 values
                       nc = numel(find(ismember(tempPairs.Y_area,'SC')));
                       outConn.SC_1(idx) = outConn.SC_1(idx) + nc/2;
                       outConn.SC_2(idx) = outConn.SC_2(idx) + nc/2;
                   end
 
               % next unit num
           end
           outConnMat.(outFieldname) = outConn;
 
       end
   end
end






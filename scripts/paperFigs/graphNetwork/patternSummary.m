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
rscTbl.X_sortVisMov(idx) = ones(numel(idx),1);
idx = find(abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Mov'},numel(idx),1);
rscTbl.X_sortVisMov(idx) = 2 + ones(numel(idx),1);
idx = find(abs(rscTbl.X_visGrade) >=2 & abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'VisMov'},numel(idx),1);
rscTbl.X_sortVisMov(idx) = 1 + ones(numel(idx),1);
% code functional types for Y
rscTbl.Y_visMovType = repmat({'Other'},size(rscTbl,1),1);
rscTbl.Y_sortVisMov = repmat(4,size(rscTbl,1),1);
idx = find(abs(rscTbl.Y_visGrade) >=2);
rscTbl.Y_visMovType(idx) = repmat({'Vis'},numel(idx),1);
rscTbl.Y_sortVisMov(idx) = ones(numel(idx),1);
idx = find(abs(rscTbl.Y_moveGrade) >=2);
rscTbl.Y_visMovType(idx) = repmat({'Mov'},numel(idx),1);
rscTbl.Y_sortVisMov(idx) = 2 + ones(numel(idx),1);
idx = find(abs(rscTbl.Y_visGrade) >=2 & abs(rscTbl.Y_moveGrade) >=2);
rscTbl.Y_visMovType(idx) = repmat({'VisMov'},numel(idx),1);
rscTbl.Y_sortVisMov(idx) = 1 + ones(numel(idx),1);
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

%% Output rsc table with only post-saccade epoch for HEB plot in R
r = rscTbl(ismember(rscTbl.epoch,'PostSaccade'),{'X_unitNum',...
                                                 'Y_unitNum',...
                                                 'X_area',...
                                                 'Y_area',...
                                                 'X_visMovType',...
                                                 'X_sortVisMov',...
                                                 'X_errGrade',...
                                                 'X_rewGrade',...
                                                 'X_poorIso',...
                                                 'Y_visMovType',...
                                                 'Y_sortVisMov',...
                                                 'Y_errGrade',...
                                                 'Y_rewGrade',...
                                                 'Y_poorIso',...
                                                 'outcome',...
                                                 'satCondition',...
                                                 'epoch',...
                                                 'rscSign',...
                                                 'signifPlusRho',...
                                                 'signifMinusRho',...
                                                 'signifRsc'...                                              
                                                 });
writetable(r,'rcode/rscTblCrossAreaPostSaccade.csv');




% Using Olaf Sporns' Brain Connectivity Toolbox for analysis
% See: https://sites.google.com/site/bctnet/
% Complex network measures of brain connectivity: Uses and interpretations.
%     Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.
%     https://doi.org/10.1016/j.neuroimage.2009.10.003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath(genpath('../BCT'))

rscTbl = getPairRscTableForHeb();
%% Get Edge features. first pair of columns are edges (from to) rest of the columns are features of Edges
allEdges = table();
allEdges.from = rscTbl.X_unitNum;
allEdges.to = rscTbl.Y_unitNum;
allEdges.fromArea = rscTbl.X_area;
allEdges.toArea = rscTbl.Y_area;
allEdges.satCondition = regexprep(rscTbl.condition,'Correct|Error.*','');
allEdges.outcome = regexprep(rscTbl.condition,'Fast|Accurate','');
allEdges.epoch = rscTbl.alignedName;
allEdges.rho = rscTbl.rhoRaw_150ms;
allEdges.pval = rscTbl.pvalRaw_150ms;


%% Get node features. Each column of this tables corresponds 
%   to the Names in the demo matrix fo the Brain Connectivity toolbox
% get List of all nodes
nodes = table();
nodes.monk = [rscTbl.monkey;rscTbl.monkey];
nodes.sessNum = [rscTbl.sessNum;rscTbl.sessNum];
nodes.sess = [rscTbl.sess;rscTbl.sess];
nodes.unitNum = [rscTbl.X_unitNum ;rscTbl.Y_unitNum];
nodes.unitArea = [rscTbl.X_area;rscTbl.Y_area];
nodes.visGrade = [rscTbl.X_visGrade;rscTbl.Y_visGrade];
nodes.visType = [rscTbl.X_visType;rscTbl.Y_visType];
nodes.moveGrade = [rscTbl.X_moveGrade;rscTbl.Y_moveGrade];
nodes.errGrade = [rscTbl.X_errGrade;rscTbl.Y_errGrade];
nodes.rewGrade = [rscTbl.X_rewGrade;rscTbl.Y_rewGrade];
nodes.poorIso = [rscTbl.X_poorIso;rscTbl.Y_poorIso];

nodes = unique(nodes,'rows');

nodes.Properties.RowNames = arrayfun(@num2str,nodes.unitNum,'UniformOutput',false);

% add additional fields to nodeFeatures (functional types)
% see rcode/plotSpkCorrHebFuncCellTypes.R #177 to #180 for this coding
nodes.sortArea = nan(size(nodes,1),1);
idx = find(ismember(nodes.unitArea,'SEF'));
nodes.sortArea(idx) = zeros(numel(idx),1);
idx = find(ismember(nodes.unitArea,'FEF'));
nodes.sortArea(idx) = zeros(numel(idx),1) + 1;
idx = find(ismember(nodes.unitArea,'SC'));
nodes.sortArea(idx) = zeros(numel(idx),1) + 2;

% see rcode/plotSpkCorrHebFuncCellTypes.R #144 to #158 for this coding
nodes.visMovType = repmat({'Other'},size(nodes,1),1);
nodes.sortVisMov = repmat(4,size(nodes,1),1);
idx = find(abs(nodes.visGrade) >= 2);
nodes.visMovType(idx) = repmat({'Vis'},numel(idx),1);
nodes.sortVisMov(idx) = ones(numel(idx),1);
idx = find(abs(nodes.moveGrade) >= 2);
nodes.visMovType(idx) = repmat({'Mov'},numel(idx),1);
nodes.sortVisMov(idx) = ones(numel(idx),1).*2;
idx = find(abs(nodes.visGrade) >= 2 & abs(nodes.moveGrade) >= 2);
nodes.visMovType(idx) = repmat({'VisMov'},numel(idx),1);
nodes.sortVisMov(idx) = ones(numel(idx),1).*3;

nodes.errorRewardType = repmat({'NA'},size(nodes,1),1);
nodes.errorRewardTxt = repmat({''},size(nodes,1),1);
idx = find(abs(nodes.errGrade) >= 2);
nodes.errorRewardType(idx) = repmat({'Error'},numel(idx),1);
nodes.errorRewardTxt(idx) = repmat({'c'},numel(idx),1);
idx = find(abs(nodes.rewGrade) >= 2);
nodes.errorRewardType(idx) = repmat({'Reward'},numel(idx),1);
nodes.errorRewardTxt(idx) = repmat({'t'},numel(idx),1);
idx = find(abs(nodes.errGrade) >= 2 & abs(nodes.rewGrade) >= 2);
nodes.errorRewardType(idx) = repmat({'Error & Reward'},numel(idx),1);
nodes.errorRewardTxt(idx) = repmat({'x'},numel(idx),1);

% Sort nodes table as needed to create the CIJ matrix
% see rcode/plotSpkCorrHebFuncCellTypes.R #181 to #184
% sortedEdges <- allEdges %>% group_by(sortArea,sortVisMov,errorRewardTxt)
% sortedEdges <- sortedEdges %>% arrange(sortArea,sortVisMov,errorRewardTxt,as.numeric(unitNum))
sortedNodes = sortrows(nodes,{'sortArea','sortVisMov','errorRewardTxt','unitNum'});


%% Build the connectivty or Adjacency matrix 
% Networks may be represented by their connectivity (adjacency) matrices.
% Rows and columns in these matrices denote nodes, while matrix entries denote links. 
CIJ = zeros(size(nodes,1));


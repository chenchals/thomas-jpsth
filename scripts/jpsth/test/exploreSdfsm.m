corrBl = load('dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrSdfs_Correct_Baseline.mat');
unitAreas = {'SEF','FEF','SC'};


%% for SEF units
allPairings = table();
unitArea = 'SEF';
currTbl = corrBl.(unitArea)(:,{'unitNum','area','condition','pairedSefUnitNums','pairedFefUnitNums','pairedScUnitNums'});
idx = find(cellfun(@(x,y,z) ~isempty(x)|~isempty(y)|~isempty(z),currTbl.pairedSefUnitNums,currTbl.pairedFefUnitNums,currTbl.pairedScUnitNums));
allPairings = [allPairings;currTbl(idx,:)]; %#ok<*FNDSB>
unitArea = 'FEF';
currTbl = corrBl.(unitArea)(:,{'unitNum','area','condition','pairedSefUnitNums','pairedFefUnitNums','pairedScUnitNums'});
idx = find(cellfun(@(x,y,z) ~isempty(x)|~isempty(y)|~isempty(z),currTbl.pairedSefUnitNums,currTbl.pairedFefUnitNums,currTbl.pairedScUnitNums));
allPairings = [allPairings;currTbl(idx,:)];
unitArea = 'SC';
currTbl = corrBl.(unitArea)(:,{'unitNum','area','condition','pairedSefUnitNums','pairedFefUnitNums','pairedScUnitNums'});
idx = find(cellfun(@(x,y,z) ~isempty(x)|~isempty(y)|~isempty(z),currTbl.pairedSefUnitNums,currTbl.pairedFefUnitNums,currTbl.pairedScUnitNums));
allPairings = [allPairings;currTbl(idx,:)];

%% Split by available parings
% all 3 pairings: SEF_SEF, SEF_FEF, SEF_SC
% unit0 area0 unit1 area1 unit2 area2 unit3 area3
%   ##    SEF   ##   SEF    ##   FEF    ##   SC
nR = size(allPairings,1);
unitPairs = table();
for ro = 1:nR
    temp = table();
    tbl = allPairings(ro,:);
    if strcmp(tbl.area,'SEF')
        unit0 = tbl.unitNum;
        area0 = tbl.area;
        % sef pairs
        t = tbl.pairedSefUnitNums;
        V =1;
        
        
    elseif strcmp(tbl.area,'FEF')
    elseif strcmp(tbl.area,'SC')
    end
end

   

%% Check o line for vis/postSacc, postRew --> is it off by 100 ms?

testSdfs = corrBl.SEF(1,:);
viz = testSdfs.Visual_sdfTsMeanStdSem{1}(:,[1,2]);
vizx = testSdfs.Visual_timeMs{1};
ps = testSdfs.PostSaccade_sdfTsMeanStdSem{1}(:,[1,2]);
psx = testSdfs.PostSaccade_timeMs{1};
pr = testSdfs.PostReward_sdfTsMeanStdSem{1}(:,[1,2]);

figure
plot(viz(:,1),viz(:,2))
hold on
plot(vizx,viz(:,2))
grid on

figure
plot(ps(:,1),ps(:,2))
hold on
plot(psx,ps(:,2))
grid on






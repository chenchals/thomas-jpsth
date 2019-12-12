function [allUnitPairs] = extractUnitsSefPairs(outcomes,epoch)
%     outcomes = {'Correct','ErrorChoice','ErrorTiming'};
%     epoch = 'Baseline';
    outFileBase = ['dataProcessed/analysis/11-18-2019/spkCorr/summary/sefPairs' epoch];
    outXlsxFile = [outFileBase '.xlsx'];    
    outMatFile = [outFileBase '.mat'];
    baseDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary';
    baseSpkCorrSdfsFile = [baseDir '/spkCorrSdfs_'];
    spkCorrRhoPvalFile = [baseDir '/spkCorrAllPairsStaticRhoPval.mat'];
    outcomeEpochs = strcat(outcomes,'_',epoch);
    SEF = 'SEF';
    FEF = 'FEF';
    SC = 'SC';
    unitAreas = {SEF,FEF,SC};

    allPairings = table();
    for oe = 1:numel(outcomeEpochs)
        outcomeEpoch = outcomeEpochs{oe};
        useCorrTbl = load([baseSpkCorrSdfsFile outcomeEpoch '.mat']);
        %% for SEF units
        for ua = 1:numel(unitAreas)
            unitArea = unitAreas{ua};
            currTbl = useCorrTbl.(unitArea)(:,{'unitNum','area','condition','pairedSefUnitNums','pairedFefUnitNums','pairedScUnitNums'});
            idx = find(cellfun(@(x,y,z) ~isempty(x)|~isempty(y)|~isempty(z),currTbl.pairedSefUnitNums,currTbl.pairedFefUnitNums,currTbl.pairedScUnitNums));
            allPairings = [allPairings;currTbl(idx,:)]; %#ok<*AGROW>
            % unitArea = 'FEF';
            % currTbl = useCorrTbl.(unitArea)(:,{'unitNum','area','condition','pairedSefUnitNums','pairedFefUnitNums','pairedScUnitNums'});
            % idx = find(cellfun(@(x,y,z) ~isempty(x)|~isempty(y)|~isempty(z),currTbl.pairedSefUnitNums,currTbl.pairedFefUnitNums,currTbl.pairedScUnitNums));
            % allPairings = [allPairings;currTbl(idx,:)];
            % unitArea = 'SC';
            % currTbl = useCorrTbl.(unitArea)(:,{'unitNum','area','condition','pairedSefUnitNums','pairedFefUnitNums','pairedScUnitNums'});
            % idx = find(cellfun(@(x,y,z) ~isempty(x)|~isempty(y)|~isempty(z),currTbl.pairedSefUnitNums,currTbl.pairedFefUnitNums,currTbl.pairedScUnitNums));
            % allPairings = [allPairings;currTbl(idx,:)];
        end
    end
    %% Split by available parings
    % all 3 pairings: SEF_SEF, SEF_FEF, SEF_SC
    % unit0 area0 unit1 area1 unit2 area2 unit3 area3
    %   ##    SEF   ##   SEF    ##   FEF    ##   SC
    nR = size(allPairings,1);
    unitPairs = table();
    for ro = 1:nR
        tbl = allPairings(ro,:);
        unit0 = 0;
        unit1 = 0;
        unit2 = 0;
        unit3 = 0;
        if strcmp(tbl.area,SEF) % possible pairs with (SEF|FEF|SC)
            % 0
            unit0 = tbl.unitNum; % only 1 SEF units possible
            % 1 paired SEF
            unit1 = tbl.pairedSefUnitNums{1};
            % 2 paired FEF
            unit2 = tbl.pairedFefUnitNums{1};
            % 3 paired SC
            unit3 = tbl.pairedScUnitNums{1};
        elseif strcmp(tbl.area,'FEF') % possible pairs with (SEF)
            % 0
            unit0 = tbl.pairedSefUnitNums{1}; % multiple SEF units possible
            % 2
            unit2 = tbl.unitNum;
        elseif strcmp(tbl.area,'SC') % possible pais with (SC)
            % 0
            unit0 = tbl.pairedSefUnitNums{1}; % multiple SEF units possible
            % 3
            unit3 = tbl.unitNum;
        end
        % pivot...
        for ii = 1:numel(unit0)
            unitPairs = [unitPairs; pivotMultipleUnits(unit0(ii),unit1,unit2,unit3)];
        end
    end
    unitPairs = unique(unitPairs,'rows');
    unitPairs.nPairedAreas = 1*(unitPairs.pairedSefUnit > 0) + 10*(unitPairs.pairedFefUnit > 0) + 100*(unitPairs.pairedScUnit > 0);
    unitPairs = sortrows(unitPairs,{'unitNum','nPairedAreas'},{'ascend','descend'});
    
    %% Add rho, pval to paired units to filter for plotting if needed
    colNames = {
        'Pair_UID'
        'X_unitNum'
        'Y_unitNum'
        'X_area'
        'Y_area'
        'XY_Dist'
        'condition'
        'alignedName'
        'rhoRaw_150ms'
        'pvalRaw_150ms'
        };
    % SEF_SEF pairs
    pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'};
    pairedUnitColNames = {'pairedSefUnit','pairedFefUnit','pairedScUnit'};
    allUnitPairs = struct();
    pairsByArea = table();
    for ii = 1:numel(pairAreas)
        pairArea = pairAreas{ii};
        pairedUnitColName = pairedUnitColNames{ii};
        paUnits = unitPairs(unitPairs.pairedSefUnit > 0,{'unitNum',pairedUnitColName});
        paUnits.Properties.VariableNames = {'unitNum','pairedUnitNum'};
        paUnits = unique(paUnits,'rows','stable');
        rhoPvals = load(spkCorrRhoPvalFile,pairArea);
        rhoPvals = rhoPvals.(pairArea);
        paUnits = innerjoin(paUnits,rhoPvals,'LeftKeys',{'unitNum','pairedUnitNum'},...
            'RightKeys',{'X_unitNum','Y_unitNum'},...
            'LeftVariables',{'unitNum','pairedUnitNum'},'RightVariables',colNames);
        pairsByArea = [pairsByArea;paUnits];
        allUnitPairs.(pairArea) = paUnits;
    end
    allUnitPairs.outcomeEpochs = outcomeEpochs;
    allUnitPairs.pairsForUnit = unitPairs;
    allUnitPairs.pairsByarea = pairsByArea;
    
    writetable(unitPairs,outXlsxFile);
    save(outMatFile,'-v7.3','-struct','allUnitPairs');
%     clearvars -except unitPairs* allPairings
end

%% Make each paired unit its own row
function [ outTbl ] = pivotMultipleUnits(unit0,u1,u2,u3)
  outTbl = table();
  if isempty(u1), u1=0; end
  if isempty(u2), u2=0; end
  if isempty(u3), u3=0; end
  
  u = [u1(:);u2(:);u3(:)];
  uAll = combnk(u,3);
  idx = ismember(uAll(:,1),u1) & ismember(uAll(:,2),u2) & ismember(uAll(:,3),u3);
  uAll = uAll(idx,:);
  uAll = sortrows(uAll,[2,3]);
  nR = size(uAll,1);
  outTbl.unitNum = repmat(unit0,nR,1);
  outTbl.unitArea = repmat('SEF',nR,1);
  outTbl.pairedSefUnit = uAll(:,1);
  outTbl.pairedFefUnit = uAll(:,2);
  outTbl.pairedScUnit = uAll(:,3);
end
   

%% Check 0 line for vis/postSacc, postRew --> is it off by 100 ms?
function [] = checkSdf(corrBl)
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
end




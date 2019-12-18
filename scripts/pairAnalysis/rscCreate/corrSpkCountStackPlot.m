%oStruct=load('temp_spkCorrStackPlotData.mat');
corrDirs = {
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-SEF/mat' 
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-FEF/mat'    
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-SC/mat'     
    'dataProcessed/analysis/spkCorr/spkCorr_FEF-FEF/mat'    
    'dataProcessed/analysis/spkCorr/spkCorr_FEF-SC/mat'     
    'dataProcessed/analysis/spkCorr/spkCorr_SC-SC/mat'  
    };
plotDirs = regexprep(corrDirs,'/mat$','/pdf');
pairAreas = regexprep(corrDirs,{'.*_','/mat','-'},{'','','_'});

for p = 1:numel(plotDirs)
    if ~exist(plotDirs{p},'dir')
        mkdir(plotDirs{p})
    end
end


%%
corrInfoFields = {
    'condition'
    'alignedName'
    'alignedEvent'
    'alignedTimeWin'
    'alignTime'
    'rasterBins'
    'rho_pval_100ms'
    };
pairInfoFields = {
    'Pair_UID'
    'X_unitNum'
    'Y_unitNum'
    'X_area'
    'Y_area'
    'isOnSameChannel'
    };
%%
if ~exist(oStruct,'var')
tic
oStruct = struct();
for d = 1:numel(corrDirs)
    %%
    pairArea = pairAreas{d};
    scFiles = dir([corrDirs{d},'/spkCorr*.mat']);
    scFiles = strcat({scFiles.folder}',filesep,{scFiles.name}');
    pdfDir = plotDirs{d};
    savePdfFlag = 1;
    %%
    %rpTemp = cell2table(cell(numel(scFiles),numel(pairInfoFields)+numel(corrInfoFields)));
    rpTemp = table();
    parfor f = 1:numel(scFiles)           
        %%
        scFile = scFiles{f};
        temp = load(scFile);
        % extract unitids for pair, rhoPval for 150ms moving win
        t = temp.spkCorr(:,corrInfoFields); 
        c = repmat(temp.cellPairInfo(1,pairInfoFields),size(t,1),1);
        rpTemp = [rpTemp; [c t]];
    end
    rpTemp.signif05 = cellfun(@(x) x(:,2)<=0.05,rpTemp.rho_pval_100ms,'UniformOutput',false);
    rpTemp.signif01 = cellfun(@(x) x(:,2)<=0.01,rpTemp.rho_pval_100ms,'UniformOutput',false);
    rpTemp.positive_05 = cellfun(@(x) x(:,1)>0 & x(:,2)<=0.05,rpTemp.rho_pval_100ms,'UniformOutput',false);
    rpTemp.negative_05 = cellfun(@(x) x(:,1)<0 & x(:,2)<=0.05,rpTemp.rho_pval_100ms,'UniformOutput',false);
    rpTemp.positive_01 = cellfun(@(x) x(:,1)>0 & x(:,2)<=0.01,rpTemp.rho_pval_100ms,'UniformOutput',false);
    rpTemp.negative_01 = cellfun(@(x) x(:,1)<0 & x(:,2)<=0.01,rpTemp.rho_pval_100ms,'UniformOutput',false);
    for r = 1:size(rpTemp,1)
        pn = zeros(size(rpTemp.positive_05{r},1),1);
        pn(rpTemp.positive_05{r}) = 1;
        pn(rpTemp.negative_05{r}) = -1;
        rpTemp.plusMinus_05{r} = pn;
        pn(:) = 0;
        pn(rpTemp.positive_01{r}) = 1;
        pn(rpTemp.negative_01{r}) = -1;
        rpTemp.plusMinus_01{r} = pn;
    end
    oStruct.(pairArea) = rpTemp;
    
end
toc
end
%%
cdataPlusMinus = [0 0 1;1 1 1;1 0 0];% idx 1=blue, 2=white, 3=red
cdataPlus = [1 1 1;1 0 0];% idx 1=white, 2=red
cdataMinus = [1 1 1;0 0 1]; % idx 1=white, 2=blue

conditionPairs = {
    {'FastErrorChoice','AccurateErrorChoice'};
    {'FastErrorTiming','AccurateErrorTiming'};
    {'FastCorrect','AccurateCorrect'};
    };
alignedNames = {
    'Baseline'
    'Visual'
    'PostSaccade'
    'PostReward'
    };

currPairArea = 'SEF_FEF';
pairedData = oStruct.(currPairArea); 

for cp = 1:numel(conditionPairs)
    condPair = conditionPairs{cp};
    for an = 1:numel(alignedNames)
        alignName = alignedNames{an};
        rowIds = find(ismember(pairedData.condition,condPair) & ismember(pairedData.alignedName,alignName));
        currData = pairedData(rowIds,:);
        currData = sortrows(currData,{'X_unitNum','Y_unitNum'});
        
    end
end







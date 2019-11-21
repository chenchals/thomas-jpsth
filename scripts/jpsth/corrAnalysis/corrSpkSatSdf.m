% plot SDFs for different epochs of units used in spike count correlation
% for paired neurons

% Use static windows
spkCorrDir = 'dataProcessed/analysis/11-18-2019/spkCorr';
fn = fullfile(spkCorrDir,'summary/spkCorrAllPairsStaticNew.mat');
%spkCorr = load(fn);
outPdfDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/sdf';
if ~exist(outPdfDir,'dir')
    mkdir(outPdfDir);
end
outPdfFn = fullfile(outPdfDir,'satSdf_');

%% Fiter data for use too plot sdfs
% Use following area pairs (these are fields of spkCorr loaded above)
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'}';
spkCorrBase = cellfun(@(x) fullfile(spkCorrDir,['spkCorr_' x],'mat'),pairAreas,'UniformOutput',false);
% Use following fields
useCols = {
    'srcFile'                      
    'pairAreas'
    'condition'
    'alignedName'
    'alignedEvent'
    'nTrials'                      
    'Pair_UID'                     
    'X_monkey'                     
    'X_sess'                       
    'X_unitNum'                    
    'Y_unitNum'                    
    'rhoRaw_200ms'                 
    'pvalRaw_200ms'                
 };
% Use significance level
useSignif = 0.05;
%%
spkCorrAllTbl = table();
for pa =1:numel(pairAreas)
    temp = spkCorr.(pairAreas{pa})(:,useCols);
    temp.srcFile = strcat(spkCorrBase{pa},filesep,temp.srcFile);
    temp = temp(temp.pvalRaw_200ms>useSignif,:);
    spkCorrAllTbl = [spkCorrAllTbl; temp];
    clear temp
end


%%
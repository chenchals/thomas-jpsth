
%% GPU or par-for CPU? GPIU is slow? parfor better
dataDir ='data'
outDir = 'data/burstAnalysis';
dataFile = 'eulsef20120910c-01.mat';
sessionNo = 14;
% Var to figure out the length of session
adVar = 'AD01';% 9622848 
% single unit variable(s) to use for burst detection
spkVar = 'DSP13b';%475415 spikes

if ~exist(outDir, 'dir')
      mkdir(outDir)
end

fName = fullfile(dataDir,dataFile);
maxTimes = whos('-file',fName,adVar);
maxTimes = maxTimes.size(1);% ms
chunkTimes = 1000*60*2; %ms
maxChunks = ceil(maxTimes/chunkTimes);
chunkAt = linspace(1,maxChunks*chunkTimes,maxChunks)';


% load spike train
spkTrain = load(fName,spkVar);
spkTrain = spkTrain.(spkVar);
tic
timeWins = arrayfun(@(x) [chunkAt(x)-1 chunkAt(x+1)],(1:numel(chunkAt)-1)','UniformOutput',false);
spkTimes = cellfun(@(x) spkTrain(spkTrain>x(1) & spkTrain<=x(2)),timeWins,'UniformOutput',false);
avgFrs = cellfun(@(x) numel(x)*1000/range(x) ,spkTimes,'UniformOutput',false);
% CPU use par-for
toc
fprintf('Progress:\n');
tic
parfor i = 1:maxChunks-1
    fprintf('\b%d, \n',i);
    oBursts{i} = poissBurst(spkTimes{i},timeWins{i}(1), timeWins{i}(2),'averageFr', avgFrs{i},'plotBursts',0);
end


outFile = fullfile(outDir,['test_1' '.mat']);
saveOutput(outFile,oBursts);

toc

so=@saveOutput;


function saveOutput(oFile, resCellArray)
   var_fx = @(fn,y) cell2mat(cellfun(@(x) x.(fn),y,'UniformOutput',false)); 
   analysisDate = datetime;
   save oFile analysisDate;
   for fn = fieldnames(resCellArray{1})
       if regexp(fn,'fieldDefinitions','matchcase')
           t.(char(fn)) = resCellArray{1}.(char(fn));           
       elseif regexp(fn,'fieldDefinitions','matchcase')
           t.(char(fn)) = resCellArray{1}.(char(fn));           
       elseif regexp(fn,'timeWin','matchcase')
           temp = var_fx(char(fn),resCellArray);
           t.(char(fn)) = [temp(1:2:end);temp(2:2:end)]';
           clear temp
       else
           t.(char(fn)) = var_fx(char(fn),resCellArray);
       end
       save oFile -append -struct t
       clearvars t
   end
end





% tic
% for i = 1:5 %maxChunks-1
%      fprintf('\b%d\n',i);
%      o3{i} = poissBurst(spkTimes{i},timeWins{i}(1), timeWins{i}(2),'plotBursts',0,'useGpu',true);
% end
% out3 = o3;
% toc
%%  quick dirty
% bobt = os1.bobT;
% eobt = os1.eobT;
% 
% ad = AD01;
% 
% 
% eegChunks = arrayfun(@(x) ad(x-20:x+20),bobt,'UniformOutput',false)';
% eegMat = mean(cell2mat(eegChunks')');
% eegMatSem = (std(cell2mat(eegChunks')'))/sqrt(size(eegChunks,1));
% plot(-20:20,eegMat)
% hold on
% plot(-20:20,eegMat-eegMatSem,'r')
% plot(-20:20,eegMat+eegMatSem,'r')
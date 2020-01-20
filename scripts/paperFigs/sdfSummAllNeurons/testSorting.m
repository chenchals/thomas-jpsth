%% Data setup
inSdfsMat = satSdfsImageTbl.Correct_PostSaccade_005.VisualSatSdfNormalized{1};
inSdfsTs = satSdfsImageTbl.Correct_PostSaccade_005.VisualTimeMs{1};

stWin = 601;
enWin = stWin + 300;
nRows = size(inSdfsMat,1);

%% Current method
[out1,outSort1] = sortSdfsMat(inSdfsMat,stWin,enWin);

imagesc(out1,[-1 1])



%% test Amir's method:
%[SigModTime, ModTime] = detectModulation(diffFunc, Mean_BaseLine, SD_BaseLine, CP1, CP2, Cdur1, Cdur2, continuityFiller);


for ii = 1:nRows
    frbl = 0;
    stdBl = 0.2;
    diffFunc = inSdfsMat(ii,:);
    CP1 = 2;
    CP2 = 3;
    Cdur1 = 50;
    Cdur2 = 50;
    continuityFiller = 20;
    
    [SigModTime, ModTime] = detectModulation(diffFunc, frBl, stdBl, CP1, CP2, Cdur1, Cdur2, continuityFiller);

end


%% Testing different sorts 
%% 1st sort by max value
absMaxIdx = arrayfun(@(x) find(abs(inSdfsMat(x,:))>=max(abs(inSdfsMat(x,:))),1),(1:nRows)');
maxVals = arrayfun(@(x) inSdfsMat(x,absMaxIdx(x)),(1:nRows)');
[maxValsOrd,maxOrd]=sort(maxVals,'descend');

%% second sort by start of vis activity, use Z
Zv = inSdfsMat(:,stWin:enWin);

[riseMs,durationMs,signOfMaxFr] = arrayfun(@(x) findRiseOf(Zv(x,:)),(1:nRows)');

[soretdVals1,sortOrd1] = sortrows([maxVals,riseMs,durationMs],[1,2,3],{'descend','ascend','descend'});

[soretdVals2,sortOrd2] = sortrows([maxVals,durationMs,riseMs],[1,2,3],{'descend','descend','ascend'});

[soretdVals3,sortOrd3] = sortrows([maxVals,riseMs],[1,2],{'descend','descend'});



% [zVTimeOrd,zVtimeOrd] = sort(startOfVis,'descend');

%%
close all
out1 = inSdfsMat(sortOrd1,:);
out2 = inSdfsMat(sortOrd2,:);
out3 = inSdfsMat(sortOrd3,:);
figure
subplot(1,4,1);imagesc(inSdfsMat(:,600:end),[-1 1]);
subplot(1,4,2);imagesc(out1(:,600:end),[-1 1]);
subplot(1,4,3);imagesc(out2(:,600:end),[-1 1]);
subplot(1,4,4);imagesc(out3(:,600:end),[-1 1]);



%% show sorted


axes(h1);imagesc(inSdfsMat,[-1 1])
axes(h2);imagesc(outSdfsSorted,[-1 1])


%% count occurances of a value on the norm. sdf
edges = -1.55:0.05:1.55;
binC = edges(1:end-1)+0.05;
stWin = 601;
enWin = 1000;

ZZ= inSdfsMat(:,stWin:enWin);
vecMin = min(ZZ(:));
vecMax = max(ZZ(:));

edges = linspace(vecMin,vecMax,100);
d = unique(diff(edges));
binC = edges(1:end-1)+d(1);


[hCounts] = arrayfun(@(x) histcounts(inSdfsMat(x,stWin:enWin),edges),(1:nRows)','UniformOutput',false);
hCounts = cell2mat(hCounts);
hCountsB = sum(hCounts.*repmat(binC,size(hCounts,1),1),2);

[riseMs,durationMs,signOfMaxFr] = arrayfun(@(x) findRiseOf(inSdfsMat(x,stWin:enWin)),(1:nRows)');


% Sortiung
[~,ordHist]=sort(hCountsB,'descend');
[~,ordRise]=sort(riseMs,'ascend');

% Now sort by risetime 
riseSdf = inSdfsMat(ordRise,:);
% Histogram sort the rise
[hCounts] = arrayfun(@(x) histcounts(riseSdf(x,stWin:enWin),edges),(1:nRows)','UniformOutput',false);
hCounts = cell2mat(hCounts);
hCountsB = sum(hCounts.*repmat(binC,size(hCounts,1),1),2);
[~,ordHist2]=sortrows([hCountsB riseMs],{'descend','ascend'});



%%
figure
subplot(1,4,1);imagesc(inSdfsMat,[-1 1]);colorbar
subplot(1,4,2);imagesc(inSdfsMat(ordHist,:),[-1 1]);colorbar
subplot(1,4,3);imagesc(inSdfsMat(ordRise,:),[-1 1]);colorbar
subplot(1,4,4);imagesc(inSdfsMat(ordHist2,:),[-1 1]);colorbar



%%
function [startTick,durationTicks,directionOfMax] = findRiseOf(vec)
   % find the sign of the max value
   maxValIdx = find(abs(vec)==max(abs(vec)));
   maxValIdx = maxValIdx(1);
   directionOfMax = 1;
   if vec(maxValIdx) < 0
       directionOfMax = -1;
   end
   d = diff(vec);
   if directionOfMax < 0
       d = diff(-vec);
   end
   startTick = find(d==max(d));
   endTick = find(d==min(d));
   startTick = startTick(1);
   endTick = endTick(end);  
   durationTicks = endTick - startTick;

end


function [ix,dypks] = findRiseAndFall(x,y)
dy = gradient(y, mean(diff(x)));                % Take Derivative
[dypks,ix] = findpeaks(dy, 'MinPeakDistance',20, 'MinPeakHeight',1E+7);
end



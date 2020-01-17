function [outSdfsSorted,outSortOrder] = sortSdfsMat(inSdfsMat,stWin,enWin)
%% sorting
%stWin = 601;
%enWin = stWin + 300;
nRows = size(inSdfsMat,1);

%% 1st sort by max value
absMaxIdx = arrayfun(@(x) find(abs(inSdfsMat(x,:))>=max(abs(inSdfsMat(x,:))),1),(1:nRows)');
maxVals = arrayfun(@(x) inSdfsMat(x,absMaxIdx(x)),(1:nRows)');
% [maxValsOrd,maxOrd]=sort(maxVals,'descend');

%% second sort by start of vis activity, use Z
Zv = inSdfsMat(:,stWin:enWin);
zVTime = arrayfun(@(x) findRiseOf(Zv(x,:))+stWin,(1:nRows)');
% [zVTimeOrd,zVtimeOrd] = sort(startOfVis,'descend');

%%
[zMaxVtm,maxTmOrd]= sortrows([maxVals zVTime],[1,2],{'descend','ascend'});

outSortOrder = maxTmOrd;
outSdfsSorted = inSdfsMat(outSortOrder,:);

end

function [idx] = findRiseOf(vec)
   d = diff(vec);
   dIdxs = arrayfun(@(x) find(d > prctile(d,x),1),(100:-1:95),'UniformOutput',false);
   dIdxs(cellfun(@isempty,dIdxs))=[];
   idx = dIdxs{1};
end
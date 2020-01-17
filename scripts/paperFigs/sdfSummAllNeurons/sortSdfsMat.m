function [outSdfsSorted,outSortOrder] = sortSdfsMat(inSdfsMat,stWin,enWin)
    % [riseMs,durationMs,signOfMaxFr] = arrayfun(@(x) findRiseOf(inSdfsMat(x,stWin:enWin)),(1:nRows)');
    outSortOrder = orderByHistogramSum(inSdfsMat,stWin,enWin);
    outSdfsSorted = inSdfsMat(outSortOrder,:);

end


function [idx] = orderByHistogramSum(inSdfsMat,stWin,enWin)
    nRows = size(inSdfsMat,1); 
    ZZ = inSdfsMat(:,stWin:enWin);
    vecMin = min(ZZ(:));
    vecMax = max(ZZ(:));

    edges = linspace(vecMin,vecMax,100);
    d = unique(diff(edges));
    binC = edges(1:end-1)+d(1);

    [hCounts] = arrayfun(@(x) histcounts(inSdfsMat(x,stWin:enWin),edges),(1:nRows)','UniformOutput',false);
    hCounts = cell2mat(hCounts);
    hCountsB = sum(hCounts.*repmat(binC,size(hCounts,1),1),2);

    [~,idx]=sort(hCountsB,'descend');
end

function [idx] = findRiseOf(vec)
   d = diff(vec);
   dIdxs = arrayfun(@(x) find(d > prctile(d,x),1),(100:-1:95),'UniformOutput',false);
   dIdxs(cellfun(@isempty,dIdxs))=[];
   idx = dIdxs{1};
end
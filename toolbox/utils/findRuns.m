function [startIdx,endIdx,runLen] = findRuns(logicalVect,zeroOrOne)
%FUNDRUNS Fing contiguous runs of 0s or 1s given a logical array
%  vect : a logical array
%  zeroOrOne : find either runs of 0 or runs of 1
% Example:
%   testVec = logical([0 0 0 0 0 1 1 1 1 0 1 0 0 1 0 1 1 1 1 0 1 0 0 0 1 1 1 0 0 0]);
%   Find runs of 1s:
%   [startIdx1,endIdx1,runLen1] = findRuns(testVec,1)
%   Find runs of 0s:
%   [startIdx0,endIdx0,runLen0] = findRuns(testVec,0)
%
    assert(islogical(logicalVect),'The input array must be logical');
    
    X = logicalVect;
    if zeroOrOne==1
        X = logicalVect - 1;
        X(X<0) = 1;
    end
    dX=diff(X);
    startIdx=find(dX==-1)+1;
    if X(1) == 0
        startIdx = [1, startIdx];
    end
    endIdx=find(dX==1);
    if X(end) == 0
        endIdx =[endIdx length(X)];
    end
    runLen = (endIdx-startIdx)+1;
end



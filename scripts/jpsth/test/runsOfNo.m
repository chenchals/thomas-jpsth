testX1 = [0 0 0 0 0 1 1 1 1 0 1 0 0 1 0 1 1 1 1 0 1 0 0 0 1 1 1 0 0 0];

X = [0 0 0 0 0 1 1 1 1 0 1 0 1 0 1 1 1 1 0 1 0 1 1 1];
%%
find0_1vec = X
[startIdx0,endIdx0,runLen0] = findRuns(find0_1vec,0)
[startIdx1,endIdx1,runLen1] = findRuns(find0_1vec,1)


%%
function [startIdx,endIdx,runLen] = findRuns(vect,trueFalse)
X = vect;
if trueFalse==1
    X = vect - 1;
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


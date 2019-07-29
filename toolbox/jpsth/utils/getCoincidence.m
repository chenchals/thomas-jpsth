function [ coinh, coinMat ] = getCoincidence(jp,lagInBins)
%GETCOINCIDENCE Summary of this function goes here
%   Detailed explanation goes here
% Note 1: Bins [1,1] is top-left and [nBins,nBins] is bottom-right
% such that the main diagonal top-left to bottom-right is synchronous
% firing of both cells at 0-lag. The length of task time is
% (nBins*binWidth).
% Note 2: Along the diagnonal is the evolution of firing (spike) synchrony
% during the task between pair of cells. Y-Cell is cell chosen for Y-Axis
% and X-Cell is cell chosen for X-Axis.
% Note 3: We can thus get evolution of firing synchrony for any lag(s), as:
% diag(0): Y-Cell fires synchronous to X-Cell with lag of (0*binWidth)
% diag(i): Y-Cell fires synchronous to X-Cell with lag of +(i*binWidth)
% diag(-i): Y-Cell fires synchronous to X-Cell with lag of -(i*binWidth)
nBins = size(jp,1);

% Aertsen
lags = -(abs(lagInBins)):abs(lagInBins);
coinMat = zeros(nBins,numel(lags));
for ii = lags   
   d = diag(jp,ii)';
   if isempty(d)
      coinMat = [];
      return
   end
   coinMat(abs(ii)+1:end,lags==ii) = d;
end  
   coinh = [(1:nBins)' sum(coinMat,2)];

end

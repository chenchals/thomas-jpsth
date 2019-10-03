% JPSTHJETCOLORMAP produce a colormap with warm colors for positive
% values and cold colors for negative values with green representing
% zero
% Modified from jpsthJetColormap.m
% jpsthJetColormap.m
% 
% Copyright 2008 Vanderbilt University.  All rights reserved.
% John Haitas, Jeremiah Cohen, and Jeff Schall

function [thisColormap] = jpsthColormap(varargin)
    
   switch numel(varargin)
       case 1
           cmLength = varargin{1};
           clim = get(gca,'CLim');
       case 2
           cmLength = varargin{1};
           clim = varargin{2};           
       otherwise
           cmLength = 65;
           clim = get(gca,'CLim');  
       
   end
	if mod(cmLength,2)==0, cmLength=cmLength+1; end
	thisColormap = zeros(cmLength, 3);
	
	cMin = clim(1);
	cMax = clim(2);
	zeroI = fix((0-cMin)/(cMax-cMin)*cmLength)+1;
	
	negLength = zeroI - 1;
	posLength = cmLength - zeroI;
	
	for i = 1:2*negLength/3
		j = i-1;
		thisColormap(i,:) = [0 j/(2*negLength/3) 1];
	end
	
	startIndex = fix(2*negLength/3)+1;
	for i = startIndex:zeroI-1
		j = length(startIndex:zeroI-1) - (i - startIndex) - 1;
		thisColormap(i,:) = [0 1 (3*j/negLength)^.7];
	end
	
	thisColormap(zeroI,:) = [0 1 0];
		 
	startIndex = zeroI+1;
	for i = startIndex:zeroI+(posLength/3)
		j = i - startIndex;
		thisColormap(i,:) = [(3*j/(posLength))^.7 1 0];
	end
	
	startIndex = fix(zeroI+(posLength/3))+1;
	for i = startIndex:cmLength
		j = length(startIndex:cmLength) - (i - startIndex) - 1;
		thisColormap(i,:) = [1 (j/(2*posLength/3)) 0];
	end
	
end

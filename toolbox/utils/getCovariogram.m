% COVARIOGRAMBRODY computes a covariogram of spike_1 and spike_2 at a 
% given lag
%
% to compile enter the following command at the MATLAB Command Window
% >> mex('covariogramBrody.c')
%
% MATLAB equivalent to covariogramBrody.c listed below
%

% covariogramBrody.m
% 
% Copyright 2008 Vanderbilt University.  All rights reserved.
% John Haitas, Jeremiah Cohen, and Jeff Schall
% 
% This file is part of JPSTH Toolbox.
% 
% JPSTH Toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% JPSTH Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with JPSTH Toolbox.  If not, see <http://www.gnu.org/licenses/>.

function [brodyCovariogram, rawCrossCorr, shuffleCorrector, sigma, sigHigh, sigLow, parts] = getCovariogram(spike1_rast, spike2_rast, psth1, psth2, psth1Sd, psth2Sd, maxLag)
    % warn user that they should use compiled version
    warning('JPSTH:compileMex', ...
    		['You are using the MATLAB version of covariogramBrody()\n' ...
				'For better performance please compile covariogramBrody.c with this command: mex(\''covariogramBrody.c\'')']);

	if nargin < 7 maxLag = 50; end
	
	trials = size(spike1_rast,1);
	trialLength =  size(spike1_rast,2);
		
	s1s2 = zeros(2*maxLag+1,1);
	p1s2 = zeros(2*maxLag+1,1);
	s1p2 = zeros(2*maxLag+1,1);
	rawCrossCorr =  zeros(2*maxLag+1,1);
	shuffleCorrector = zeros(2*maxLag+1,1);
	%ij = zeros(4,2*maxLag+1);
	for ii = 1:2*maxLag+1
		currentLag = ii - maxLag - 1;
		
		if currentLag < 0 
			jVector = 1-currentLag:trialLength;
		else 
			jVector = 1:trialLength-currentLag; 
		end
		
		for jj = jVector
			rawCrossCorr(ii) = rawCrossCorr(ii) + mean(spike1_rast(:,jj) .* spike2_rast(:,jj+currentLag));
			shuffleCorrector(ii) = shuffleCorrector(ii) + psth1(jj) * psth2(jj+currentLag);
			s1s2(ii) = s1s2(ii) + psth1Sd(jj).^2 * psth2Sd(jj+currentLag).^2;
			p1s2(ii) = p1s2(ii) + psth1(jj).^2 * psth2Sd(jj+currentLag).^2;
			s1p2(ii) = s1p2(ii) + psth1Sd(jj).^2 * psth2(jj+currentLag).^2;
            %ij(:,i) =[currentLag;i;j;j+currentLag];
		end
	end
	
	brodyCovariogram = rawCrossCorr - shuffleCorrector;
	
	sigma = sqrt((s1s2 + p1s2 + s1p2)/trials);
	sigHigh = 2 * sigma;
	sigLow = -2 * sigma;
    
    parts.rawCrossCorr = rawCrossCorr;
    parts.shuffleCorrector = shuffleCorrector;
    parts.s1s2 = s1s2;
    parts.p1s2 = p1s2;
    parts.s1p2 = s1p2;
    parts.sigma = sigma;
    %parts.ij = ij;
end

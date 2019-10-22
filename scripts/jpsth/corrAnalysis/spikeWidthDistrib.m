% Plot spike width distribution for 
% SEF (D,E), FEF (D), and SC (D,E)

% get ninfo for unitNum corresponding to units of SEF, FEF, SC from Da and
% Eu. Use ninfo_nstats_SAT.mat
ninfoFile = 'dataProcessed/dataset/ninfo_nstats_SAT.mat';
ninfo = load(ninfoFile);
ninfo = struct2table(ninfo.ninfo);
wavWidthTbl = ninfo(~cellfun(@isempty, regexp(ninfo.monkey,'D|E','match')) ...
                  & ~cellfun(@isempty, regexp(ninfo.area,'^(SEF|FEF|SC)$','match')),...
                    {'monkey','unitNum','area'});
wavWidthTbl.unitNum = double(wavWidthTbl.unitNum);

% get spike widths for search from Unit_ddd.mat files for unitNums in
% wavWidthTbl. The colName to get is 'wavWidthSearch'.  The width is in
% sample space where each unit corresponds to 1/40000 sec or 25 microsecs
sampTime = 0.025; %ms
wavDir = 'dataProcessed/dataset/wavesNew';

nUnits = size(wavWidthTbl,1);
wavWidths = struct();
parfor jj = 1:nUnits
    unitNum = wavWidthTbl.unitNum(jj);
    unitStr = num2str(unitNum,'Unit_%03d');
    temp = load(fullfile(wavDir,[unitStr '.mat']));
    wavWidths(jj).unitNum = unitNum;
    wavWidths(jj).wavWidthSearch = temp.(unitStr).wavWidthSearch;
end

wavWidthsTbl = innerjoin(wavWidthTbl, struct2table(wavWidths));

wavWidthsTbl.wavWidthSearchMs = cellfun(@(x) sampTime*abs(cell2mat(x)), wavWidthsTbl.wavWidthSearch,'UniformOutput',false);

wavWidthsTbl.wavWidthSearchMean = cellfun(@mean, wavWidthsTbl.wavWidthSearchMs);
wavWidthsTbl.wavWidthSearchStd = cellfun(@std, wavWidthsTbl.wavWidthSearchMs);

              
wavWidthsTblStats = grpstats(wavWidthsTbl(:,{'monkey','area','wavWidthSearchMean'}),...
                  {'monkey','area'},{'mean','std'},'DataVars',{'wavWidthSearchMean'});
              
sefWidths = cell2mat(wavWidthsTbl.wavWidthSearchMs(strcmp(wavWidthsTbl.area,'SEF')));
fefWidths = cell2mat(wavWidthsTbl.wavWidthSearchMs(strcmp(wavWidthsTbl.area,'FEF')));
scWidths = cell2mat(wavWidthsTbl.wavWidthSearchMs(strcmp(wavWidthsTbl.area,'SC')));

minWidth = min([sefWidths;fefWidths;scWidths]);
maxWidth = max([sefWidths;fefWidths;scWidths]);


binLims = [0 max([max(xWavWidths),max(yWavWidths)])+100];





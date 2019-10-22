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

binWidth = 0.025;
[sefC,be] = histcounts(sefWidths,'BinLimits',[0 1],'BinWidth',binWidth);
[fefC] = histcounts(fefWidths,'BinLimits',[0 1],'BinWidth',binWidth);
[scC] = histcounts(scWidths,'BinLimits',[0 1],'BinWidth',binWidth);

xLims = [0.05 0.7];
units = wavWidthsTblStats.GroupCount;
means = round(wavWidthsTblStats.mean_wavWidthSearchMean.*1000);
stds = round(wavWidthsTblStats.std_wavWidthSearchMean.*1000);

h_fig = figure;
%% SEF
subplot(3,1,1)
h_sef = bar(be(2:end),sefC./sum(sefC),'BarWidth',1,'EdgeColor','none','FaceAlpha',0.5);
xlim(xLims)
xlabel('Spike width (ms)')
ylabel('Normalized counts')
title('SEF Spike width distribution')
idxD = find(strcmp(wavWidthsTblStats.monkey,'D') & strcmp(wavWidthsTblStats.area,'SEF'));
idxE = find(strcmp(wavWidthsTblStats.monkey,'E') & strcmp(wavWidthsTblStats.area,'SEF'));
txt = {sprintf('%s(n=%d) %3d\\pm%2d', 'D', units(idxD), means(idxD), stds(idxD));...
       sprintf('%s(n=%d) %3d\\pm%2d', 'E', units(idxE), means(idxE), stds(idxE))};
yl=get(gca,'YLim');
h_tsef =text(min(xLims)+0.05*range(xLims),max(yl)-0.25*range(yl),txt,'Interpreter','tex','FontWeight','bold');
%% FEF only Darwin
subplot(3,1,2)
h_fef = bar(be(2:end),fefC./sum(fefC),'BarWidth',1,'EdgeColor','none','FaceAlpha',0.5);
xlim(xLims)
xlabel('Spike width (ms)')
ylabel('Normalized counts')
title('FEF Spike width distribution')
idxD = find(strcmp(wavWidthsTblStats.monkey,'D') & strcmp(wavWidthsTblStats.area,'FEF'));
txt = {sprintf('%s(n=%d) %3d\\pm%2d', 'D', units(idxD), means(idxD), stds(idxD));...
       'E(n=0)'};
yl=get(gca,'YLim');
h_tfef =text(min(xLims)+0.05*range(xLims),max(yl)-0.25*range(yl),txt,'Interpreter','tex','FontWeight','bold');

%% SC
subplot(3,1,3)
h_sc = bar(be(2:end),scC./sum(scC),'BarWidth',1,'EdgeColor','none','FaceAlpha',0.5);
xlim(xLims)
xlabel('Spike width (ms)')
ylabel('Normalized counts')
title('SC Spike width distribution')
idxD = find(strcmp(wavWidthsTblStats.monkey,'D') & strcmp(wavWidthsTblStats.area,'SC'));
idxE = find(strcmp(wavWidthsTblStats.monkey,'E') & strcmp(wavWidthsTblStats.area,'SC'));
txt = {sprintf('%s(n=%d) %3d\\pm%2d', 'D', units(idxD), means(idxD), stds(idxD));...
       sprintf('%s(n=%d) %3d\\pm%2d', 'E', units(idxE), means(idxE), stds(idxE))};
yl=get(gca,'YLim');
h_tsc =text(min(xLims)+0.05*range(xLims),max(yl)-0.25*range(yl),txt,'Interpreter','tex','FontWeight','bold');






function [ sessionFiles ] = createSatMat2PlxMap()
    % create a mapping file from mat file to plx file for SAT: note only
    % -RH_SEARCH.mat files are used
    monkIds = containers.Map({'D','E','Q','S'},{'Darwin','Euler','Quincy','Seymour'});
    % load ninfo
    ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat');
    ninfo = struct2table(ninfo.ninfo);
    sessionFiles = unique(ninfo(:,{'sess','sessNum','monkey'}),'rows');
    % for Quincy and Seymour, the plx filename pattern = [Q|S]mmddyy00[1-3]
    fx_morphDate = @(x) regexprep(x,'([Q|S])\d{2}(\d{2})(\d{2})(\d{2})(\d{3})','$1$3$4$2$5');
    sessionFiles.plxSessName2Use = cellfun(fx_morphDate,sessionFiles.sess,'UniformOutput',false);
    % assume mat datafile name and check if exists
    sessionFiles.matDatafile = arrayfun(@(x) regexprep(fullfile('data',monkIds(sessionFiles.monkey{x}),...
                             'SAT/Matlab',[sessionFiles.sess{x} '-RH_SEARCH.mat']),'\','/'),...
                             (1:size(sessionFiles,1))','UniformOutput',false);
    sessionFiles.existMatDatafile = cellfun(@(x) exist(x,'file')==2,sessionFiles.matDatafile);
    % assume plx datafile name (in sorted dir) and check if exists
    sessionFiles.plxDatafile = arrayfun(@(x) regexprep(fullfile('data',monkIds(sessionFiles.monkey{x}),...
                             'SAT/Plexon/Sorted',[sessionFiles.plxSessName2Use{x} '-RH.plx']),'\','/'),...
                             (1:size(sessionFiles,1))','UniformOutput',false);
    sessionFiles.existPlxDatafile = cellfun(@(x) exist(x,'file')==2,sessionFiles.plxDatafile);
    % fix non-existent plexon files.... check if a similar file exist in
    % plexon/Unsorted/xxxx.plx (no '-RH' extention)
    plxNotExistIds= find(sessionFiles.existPlxDatafile==0);
    for ii = plxNotExistIds
        fn = sessionFiles.plxDatafile{ii};
        fnNew = regexprep(fn,'-RH','');
        if exist(fnNew,'file')
            sessionFiles.plxDatafile(ii) = fnNew;
        else 
            fnNew = regexprep(fnNew,'Sorted','Unsorted');
            if exist(fnNew,'file')
                sessionFiles.plxDatafile(ii) = {fnNew};
            end
        end
    end
    sessionFiles.existPlxDatafile = cellfun(@(x) exist(x,'file')==2,sessionFiles.plxDatafile);
    save('dataProcessed/dataset/SessionFiles_SAT.mat','sessionFiles');
end
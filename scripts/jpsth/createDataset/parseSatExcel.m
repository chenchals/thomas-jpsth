function [excelInfos] = parseSatExcel(excelSummaryFile)
% PARSESATEXCEL is specific for processing xlsx files that are recoded
% versions of Rich's xlsx files that cummarize SAT data in 
% D, E, Q, and S.
% Note session name for Q and S is different for plx and.mat files
% PLX file: Qmmddyy00n.plx --> Q20yymmdd00n-RH*.mat
% 
    [~,~,rawCell] = xlsread(excelSummaryFile);
    % find numbers
    containsNumbers = cellfun(@isnumeric,rawCell);
    %# convert to string
    rawCell(containsNumbers) = cellfun(@num2str,rawCell(containsNumbers),'UniformOutput',false);
    % Single unit quality metric
    % Not used...
    temp = rawCell(1:4,1); % top 4 rows
    qualityMetric = struct();
    for ii = 1:size(temp,1)
        raw = temp{ii};
        raw_ = split(regexprep(raw,'\s*=\s*','='),'=');
        qualityMetric(ii).raw = char(raw);
        qualityMetric(ii).val = str2double(raw_{1});
        qualityMetric(ii).comment = char(raw_{2});
        if isnan(qualityMetric(ii).val)
            % convert '< .5', was always 0.25
            qualityMetric(ii).val = 0.25;
        end
    end
    clearvars ii temp raw raw_;
    lastColumn = 23;
    colNames = rawCell(6,1:lastColumn);
    temp = rawCell(7:end,1:lastColumn);
    allNanRows = cell2mat(arrayfun(@(x) all(isnan([temp{x,:}])),(1:size(temp))','UniformOutput',false));
    temp(allNanRows,:)=[];
    sessionNames = regexp([temp{:,1}],'[A-Z]\d+','match')';
    expr = char(join(sessionNames,'|'));
    sessionNameRows = find(cell2mat(arrayfun(@(x) max([regexp(char(temp{x,1}),expr,'start'),0]) == 1,...
        (1:size(temp,1))','UniformOutput',false)));
    nSessions = numel(sessionNameRows);
    nRows = size(temp,1);
    excelInfos = [];
    for ii = 1:nSessions
        if ii < nSessions
            sessionRows = sessionNameRows(ii):sessionNameRows(ii+1)-2;
        else
            sessionRows = sessionNameRows(ii):nRows;
        end
        % Process session Information
        tempSession = cell2table(temp(sessionRows(1),1:2),'VariableNames',{'SessionName', 'SessionNotes'});
        tempTaskTrials = cell2table(temp(sessionRows(2:4),1:2),'VariableNames',{'TaskName', 'NTrials'});
        tempTaskTrials = cell2table(tempTaskTrials{:,2}','VariableNames',strcat('nTrials_',tempTaskTrials{:,1})');
        % process Unit information
        tempUnitTable = cell2table(temp(sessionRows(2:end),3:end),'VariableNames',colNames(3:end));
        tempUnitTable(strcmpi(tempUnitTable.Unit,'NaN'),:) = [];
        tempUnitTable.UnitFormatted = regexprep(tempUnitTable.Unit,'^(\d[a-z])','0$1');
        tempUnitTable.UnitName = strcat('DSP',tempUnitTable.UnitFormatted);
        %% The plx sessions were renamed after translation for: Q and S
        % Qmmddyy00n.plx --> Q20yymmdd00n.mat
        tempSession.MatSessionName = tempSession.SessionName;
        for c = 1: size(tempSession,1)
            plxName = tempSession.MatSessionName{c};
            if startsWith(plxName,'Q') || startsWith(plxName,'S')
                matName = [plxName(1) '20' plxName(6:7) plxName(2:5) plxName(8:end)];
                tempSession.MatSessionName{c} = matName;
            end
        end
        
        for c = 1:size(tempUnitTable,1)
            excelInfos = [excelInfos; [tempSession tempTaskTrials tempUnitTable(c,:)]]; %#ok<AGROW>
        end
    end
end


function plotAddPairInfo(H_axes, unitInfosTable)
%ADDCELLPAIRINFO Summary of this function goes here
%   Detailed explanation goes here

    fs = 8;
    if ismac
        fs = 10;
    end

    axes(H_axes)
    infoTable = unitInfosTable;
    varNames=infoTable.Properties.VariableNames;
    
    nRows = size(infoTable,1);
    
    % if it is a row from JpsthPairCellInfoDB
    if contains(varNames,'Pair_UID')
        sessionNotes = infoTable.SessionNotes;
        infoTable.MatSessionName = [];
        infoTable.SessionNotes = [];
    end
    
    varNames(contains(varNames,'trRem'))=[];
    fontSize = fs;
    columnGap = 0.005;
    nCols = numel(varNames);
    % we are plotting this from bottom
    xPos=0.005;
    yHPos = 0.98;
    for c = 1:nCols
        colName = varNames{c};
        hText{1} = text(xPos,yHPos,colName,'Interpreter','none',...
            'FontWeight','bold','FontSize',fontSize,...
            'VerticalAlignment', 'top','HorizontalAlignment','left');
        vals = infoTable.(colName);
        valTxt = cell(2,1);
        for ii = 1:nRows
            if iscell(vals)
                value = vals{ii};
            else
                value = vals(ii);
                if islogical(value)
                    value = double(value);
                end
            end
            if isnumeric(value)             
                 if numel(value) > 1
                     valTxt{ii,1} = ['[' num2str(value) ']'];
                 else
                     valTxt{ii,1} = num2str(value);
                 end

            else
                valTxt{ii,1} = char(value);
            end
        end
        hText{2} = text(xPos,yHPos-0.26,valTxt,'Interpreter','none',...
            'FontWeight','normal','FontSize',fontSize,...
            'VerticalAlignment', 'top','HorizontalAlignment','left');
        xPos = getNextXPos(hText, columnGap);
    end
    % now for session notes 
    if contains(varNames,'SessionNotes')
         addSessionNotes(xPos,yHPos,fontSize,sessionNotes);
    end
    
    % write Analysis date
    annotation('textbox','Position',[0.01 0.012 0.16 0.01],'String', char(datetime),...
        'FontSize',fontSize*0.75,'FitBoxToText','off','EdgeColor','none','Interpreter','none');

end

%% Get next X position fron the previous plot extens
function [ xPos ] = getNextXPos(hText, columnGap)
    xtnt = cell2mat(cellfun(@(x) get(x,'Extent'),hText','UniformOutput',false));
    xPos = max(xtnt(:,1)) + max(xtnt(:,3)) + columnGap;
end

%% add session notes
function [] = addSessionNotes(xPos,yHPos,fontSize,sessionNotes)
        text(xPos,yHPos,'SessionNotes','Interpreter','none',...
            'FontWeight','bold','FontSize',fontSize,...
            'VerticalAlignment', 'top','HorizontalAlignment','left');
        if numel(unique(sessionNotes)) > 1
            fontSize = fontSize*0.6;
        end
        annotation('textbox','Position',[xPos-0.01 0.0125 0.16 0.045],'String', unique(sessionNotes),...
            'FontSize',fontSize,'FitBoxToText','off','EdgeColor','none');


end


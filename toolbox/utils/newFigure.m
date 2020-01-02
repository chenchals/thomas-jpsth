function [H_Figure] = newFigure()
%NEWFIGURE Create a new figure with screen coordinates
% You can now use saveFigPdf to correctly save the figure
% see also SAVEFIGPDF
%
    fs = 6;
    if ismac
        fs = 8;
    end
    set(0,'units','pixels');
    set(0,'defaulttextfontsize',fs,...
        'defaulttextfontname','Arial',...
        'defaultaxesfontsize',fs,...
        'defaultaxeslinewidth',0.05);
    margin = 10; %pixels

    ss = [20 20 1500 990];
    FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
    %Main figure window
    H_Figure=figure('Position',FigPos,...
        'color',[1 1 1],'numbertitle','off','renderer','opengl',...
        'renderermode','manual','menubar','none',...
        'Tag','H_SaveFigPdf');
    orient landscape
    set(H_Figure,'units','normalized')
end
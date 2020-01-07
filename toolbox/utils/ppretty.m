function [] = ppretty(image_size, varargin)
% PPRETTY `prettify' a figure
%   [] = PPRETY(varargin) beautifies a figure for use in publications.
%   A number of figure properties are changed. Specifically the size of the
%   figure is change to fit in the proper orientation on a US-legal sized
%   piece of paper (the figure is centered inside the page). In addition, 
%   the fonts on each of the axes are changed to Arial size 16. The
%   background of the figure is also set to white (rather than the default
%   grey).
%
%   PPRETTY takes a number of optional arguments including
%     1) `margin': the margin size in inches
%     2) `square': A flag that resizes the paper to make it a 10x10 square.
%     3) `cut': A flag that says whether shading should be cut at the edge
%     of the figure
%
args = getopt(varargin, {'yRight',{'font=',8}, {'margin=', 0.5}});

% Default values
set(gcf, 'DefaultTextFontSize', args.font); % [pt]
set(gcf, 'DefaultAxesFontSize', args.font); % [pt]
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');

% Figure Values
set(gcf, 'Color', [1.0 1.0 1.0]);
set(gcf, 'Units', 'inches');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', image_size + [1 1] * 2 * args.margin);
set(gcf, 'PaperOrientation', 'portrait');

% Set the position of the figure
paper_size = get(gcf, 'PaperSize');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'Clipping', 'off');
set(gcf, 'PaperPosition', [args.margin, args.margin, paper_size(1), paper_size(2)]);
pos = get(gcf, 'Position');
% pos(1:2) = pos(1:2) + [15,-4];
pos(3:4) = [paper_size(1) paper_size(2)];
set(gcf, 'Position', pos);

% Make sure we use correct renderer
set(gcf, 'Renderer', 'painters');

% Axis Values
all_axes = findall(gcf,'type','axes');
for i = 1:length(all_axes)
    curr_axis = all_axes(i);
    set(curr_axis, 'FontName', 'Arial' );
    set(curr_axis, 'FontUnits', 'points');
    set(curr_axis, 'Box', 'off');
    set(curr_axis, 'TickDir', 'out');
    set(curr_axis, 'TickLength', [0.02 0.02]);
    set(curr_axis, 'XMinorTick', 'on');
    set(curr_axis, 'YMinorTick', 'on');
%     set(curr_axis, 'XColor', [0 0 0]);
%     set(curr_axis, 'YColor', [0 0 0]);
    set(curr_axis, 'LineWidth', 0.75);
    if (args.yRight)
      set(curr_axis, 'YAxisLocation','right')
    end
end

% Change all text sizes
all_text = findall(gcf, 'type', 'text');
for i = 1:length(all_text)
    set(all_text(i), 'FontName', 'Arial');
end

% Change the legends so that there aren't any boxes
all_legends = findall(gcf, 'tag', 'legend');
for i = 1:length(all_legends)
    set(all_legends(i), 'box', 'off');
end

end

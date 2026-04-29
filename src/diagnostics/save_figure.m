function save_figure(fig, filepath_without_ext, dpi)
%SAVE_FIGURE Saves a figure as both high-resolution PNG and editable FIG.
%
%   Args:
%       fig (figure handle): The figure to save.
%       filepath_without_ext (char): Full path without file extension.
%       dpi (int, optional): Resolution in dots per inch. Default 300.

    if nargin < 3
        dpi = 300;
    end

    % Move all legends to best position outside plot if possible
    legends = findall(fig, 'Type', 'Legend');
    for i = 1:numel(legends)
        legends(i).Location = 'eastoutside';
        legends(i).FontSize = 7;
    end

    % Save high-resolution PNG
    print(fig, sprintf('%s.png', filepath_without_ext), '-dpng', sprintf('-r%d', dpi));

    % Save editable FIG
    savefig(fig, sprintf('%s.fig', filepath_without_ext));
end
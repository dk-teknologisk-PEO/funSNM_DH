function save_compressed_png(fig, filepath, dpi)
%SAVE_COMPRESSED_PNG Saves a figure as a smaller PNG file.
%   Uses lower DPI and JPEG-like compression via print().
%
%   Args:
%       fig: figure handle
%       filepath: output path (e.g., 'results/plot.png')
%       dpi: resolution (default 100; use 72 for smallest, 150 for decent quality)

    if nargin < 3
        dpi = 100; % down from MATLAB's default 150-300
    end
    
    % Reduce figure size before saving
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'Position', [100 100 900 600]); % smaller canvas (pixels)
    
    print(fig, filepath, '-dpng', sprintf('-r%d', dpi));
end
function my_export_fig(varargin)
% A wrapper for 'export_fig' that ensures that the white background in enabled. 

set(gcf, 'Color', 'w');
export_fig(varargin{:});

end


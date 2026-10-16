function folder=figure_folder(root,varargin)
% Use the split-format figure layout when present; otherwise use the workspace.
if isfolder(fullfile(root,'figures_matlab')),base='figures_matlab';else,base='figures';end
folder=fullfile(root,base,varargin{:});
end

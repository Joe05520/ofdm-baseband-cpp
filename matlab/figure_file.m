function path=figure_file(folder,varargin)
% Resolve PNG, SVG and FIG files in either supported folder layout.
path=char(fullfile(folder,varargin{:}));marker=[filesep 'figures_matlab' filesep];at=strfind(path,marker);
if ~isempty(at)
    at=at(end);base=path(1:at-1);relative=path(at+length(marker):end);
    [~,~,extension]=fileparts(relative);
    path=fullfile(base,'figures_matlab',extension(2:end),relative);
end
parent=fileparts(path);if ~isfolder(parent),mkdir(parent);end
end

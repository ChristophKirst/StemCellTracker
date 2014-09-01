function s = cellfunc(varargin)
%
% s = cellfunc(varargin)
%
% description:
%    wroks as cell fun but by default returns cell
%

s = cellfun(varargin{:}, 'UniformOutput', false);

end
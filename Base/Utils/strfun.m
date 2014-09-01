function s = strfun(varargin)
%
% s = strfun(varargin)
%
% description:
%    wroks as cell fun but by default returns cell
%

s = cellfun(varargin{:}, 'UniformOutput', false);

end
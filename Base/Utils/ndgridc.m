function grc = ndgridc(varargin)
%
% varargout = ndgrid(varargin)
%
% description:
%    ndgrid for cells but returning a cell containing the individual coordinate values of the grid
%
% See also: ndgrid

if nargin==0
   error('ndgridc: not enough input arguments');    
end

if nargin==1
    grc = varargin;
    return
end


j = 1:nargin;
siz = cellfun(@numel,varargin);

grc = cell(1,nargin);

for i=1:nargin
   x = varargin{i};
   s = ones(1,nargin);
   s(i) = numel(x);
   x = reshape(x,s);
   s = siz;
   s(i) = 1;
   grc{i} = repmat(x,s);
end


end
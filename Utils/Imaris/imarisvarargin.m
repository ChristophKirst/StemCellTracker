function [imaris, var, nin] = imarisvarargin(varargin)
%
% [imaris, var, nin] = imarisvarargin(varargin)
%
% description:
%    checks if first argument is imaris instance and returns this and rest of input
%    
% output:
%    imaris    reference to imaris application
%    var       remainder of varargin
%    nin       length of remainder of varargin
%
% See also: imarisinstance

if nargin < 1
   imaris = imarisinstance();
   var = {};
elseif isimarisid(varargin{1})
   imaris = imarisinstance(varargin{1});
   var = varargin(2:end);
else
   imaris = imarisinstance();
   var = {varargin{:}}; %#ok<CCAT1>
end

if nargout > 2
    nin = length(var);
end
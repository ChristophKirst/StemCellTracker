function [imaris, var] = imarisvarargin(varargin)
%
% [imaris, var] = imarisvarargin(varargin)
%
% description:
%    checks if first argument is imaris instance and returns this and rest of input
%    

if nargin < 1
   imaris = imarisinstance();
   var = {};
elseif isimarisid(varargin{1})
   imaris = imarisinstance(varargin{1});
   var = varargin(2:end);
else
   imaris = imarisinstance();
   var = {varargin};
end
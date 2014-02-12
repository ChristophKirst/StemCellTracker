function [imaris, varargout] = imarisvarargin(varargin)
%
% [imaris, varargout] = imarisvarargin(varargin)
%
% description:
%    checks if first argument is imaris instance and returns this and rest of input
%    

varargout = {};

if nargin < 1
   imaris = imarisinstance();
elseif isimarisid(varargin{1})
   imaris = imarisinstance(varargin{1});
   varargout = {varargin{2:end}};
else
   imaris = imarisinstance();
end
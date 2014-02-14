function [imaris, varout, nin] = imarisvarargin(varin)
%
% [imaris, var, nin] = imarisvarargin(varin)
%
% description:
%    checks if first argument is imaris instance and returns this and rest of input
%
% input
%    varin     varargin of top function
%
% output:
%    imaris    reference to imaris application
%    varout    remainder of varargin
%    nin       length of remainder of varargin
%
% note: varin is a cell of inputs of length 1
%
% See also: imarisinstance

nargin = length(varin);

if nargin < 1
   imaris = imarisinstance();
   varout = {};
elseif isimaris(varin{1})
   imaris = varin{1};
   %varout = { varin{2:end} };
   varout = varin(2:end);
elseif isimarishead(varin{1})
   imaris = imarisinstance(varin{1});
   %varout = { varin{2:end} };
   varout = varin(2:end);  
else
   imaris = imarisinstance();
   %varout = {varin{:}};
   varout = varin;
end

if nargout > 2
    nin = length(varout);
end
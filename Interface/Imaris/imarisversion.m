function version = imarisversion(varargin)
%
% version = imarisversion(imaris)
%
% description:
%    returns version of running imaris
%
% input:
%    imaris   (optional) Imaris Application instance
%
% output:
%    version   Imaris version
%
% See also: imarisstart

imaris = imarisvarargin(varargin);
version = imaris.GetVersion;

end





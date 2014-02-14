function version = imarisversion(varargin)
%
% version = imarisversion(imaris)
%
% description:
%    returns version of running imaris
%
% input:
%    imaris   (optional) Imaris Application instance

imaris = imarisvarargin(varargin);
version = imaris.GetVersion;

end





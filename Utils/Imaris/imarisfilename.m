function filename = imarisfilename(varargin)
%
%  filename = imarisfilename(imaris)
%
% description:
%    returns current open file name in imaris
%
% input:
%    imaris   (optional) Imaris Application instance

imaris = imarisvarargin(varargin{:});
filename = imaris.GetCurrentFileName;

end


  
  
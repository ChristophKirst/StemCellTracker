function img = imreadStack(filename, varargin)
%
% stack = imreadStack(filename, varargin)
%
% description:
%     reads a image as a stack (tpycially tiff)
%
% input:
%     filename   file name of the image
%     varargin   inputs to imreadBF
%
% output:
%     img        the image stack
%
% See also: imreadBF

img = imreadBF(filename, varargin{:});
%img = imfrmtReformat(img, 'XYZC', 'YXZC');

end
   
function imwriteBF(img, varargin)
%
% for consistence of input and output formats and plotting
%

img = imfrmtReformat(img, 'XYZLCT' , 'yXZCT');
imwrite(img, varargin{:})

end
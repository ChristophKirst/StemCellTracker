function imwrite_bf(img, varargin)
%
% for consistence of input and output formats and plotting
%

img = impqlpermute(img, 'pqlct' , 'yplct');
imwrite(img, varargin{:})

end
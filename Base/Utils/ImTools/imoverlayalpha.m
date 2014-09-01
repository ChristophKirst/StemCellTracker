function imgo = imoverlayalpha(img, imgi, varargin)
%
% imgo = imoverlayalpha(img, imgo)
% imgo = imoverlayalpha(img, imgo, alpha)
%
% description:
%    overlays image img with imgo using the transparency values alpha
%
% input:
%    img    rgb image
%    imgi   rgb image to overlay
%    alpha  (optional) alpha channel (alpha = max(imgo,3))
%
% output:
%    imgo   overlayed gb image

if nargin < 3
   alpha = max(imgi,[], 3);
   alpha = cat(3, alpha, alpha, alpha);
else
   if isscalar(varargin{1})
      alpha = max(imgi,[], 3); 
      alpha = mat2gray(alpha);
      alpha = varargin{1} * cat(3, alpha, alpha, alpha);
   else
      alpha = varargin{:};
      alpha = cat(3, alpha, alpha, alpha);
   end
end

if ~isequal(size(imgi), size(alpha)) || ~isequal(size(img), size(alpha))
   size(imgi)
   size(alpha)
   size(img)
   
   error('imoverlayalpha: inconsistent image and alpha channel dimensions');
end

imgo = img .* (1 - alpha) + imgi .* alpha;

end

function imgr = imresample(img, fac)
%
% stack = imresample(img, fac)
%
% description:
%    resmaples the image img
%
% input:
%    img        input image
%    fac        resampling factor  number of [facX, facY (, facZ)]
%
% output:
%    imgr       resmapled image
%

% ceck for color
imf = imfrmtFormat(img);
dim = length(imf);
cdim = strfind(imf, 'c');
if ~isempty(cdim)
   cdim = cdim(1);
   sdim = dim - 1;
else
   sdim = dim;
end

if isscalar(fac)
   fac = ones(1,dim) * fac;
else
   if length(fac) ~= sdim
      error('imresample: image dimension and resampling factor do not agree');
   end
end

resamp = cell(1,dim);
isiz = size(img);

p = 1;
for i = 1:dim
   if ~isempty(cdim) && i == cdim
      resamp{i} = 1:isiz(i);
   else
      resamp{i} = 1:fac(p):isiz(i);
      p = p + 1;
   end
end

imgr = img(resamp{:});

end



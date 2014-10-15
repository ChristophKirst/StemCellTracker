function [bkg, flt] = backgroundAndFlatfieldFromMinAndMean(imgs, varargin)
%
% [bkg, flt] = backgroundAndFlatfieldFromMinAndMean(imgs, varargin)
%
% description:
%    take minimum over images imgs to estimate backgournd
%
% input:
%    imgs    cell array of images, or ImageSource.celldata
%    param   (optional) parameter struct with entries
%            .filter          filter estimated background with Gaussian (true)
%            .filtersize      Gaussian filter size after taking the minimum (=100)
%            .maximages       maximal number of images (=500)
%            range specs
%
% output:
%    bkg     estimated background image 
%
% See also: backgroundFromOpening

param = parseParameter(varargin);

nmax = getParameter(param, 'maximages', 500);

if iscell(imgs)
   imgs = imgs(:);
   n = length(imgs);
   if n > nmax
      ids = randperm(n, nmax);
      n = nmax;
   else
      ids = 1:n;
   end
   isIS = false;
elseif isa(imgs, 'ImageSource')
   ids = imgs.cellIndex(param);
   ids = ids(:);
   n = length(ids);
   if n > nmax
      ids = ids(randperm(n, nmax));
      n = namx;
   end
   isIS = true;
end

if n == 0
   error('backgroundAndFlatfieldFromMinAndMean: not images to estimate background from!');
end
   
if nargout > 1
   flatfield = true;
else
   flatfield = flase;
end

for i = 1:n
   if isIS
      img = imgs.data(ids(i));
   else
      img = imgs{ids(i)};
   end

   if i == 1
      bkg = img;
      if flatfield
         flt = img;
      end
   else
      bkg = min(bkg, img);
      if flatfield
         flt = ((i-1)*flt + img) / i;
      end
   end
end

if getParameter(param, 'filter', true)
   fs = getParameter(param, 'filtersize', 100);
   bkg = filterGaussian(bkg, fs);
   if flatfield
      flt = filterGaussian(bkg,fs);
   end
end

end
function bkg = backgroundFromMin(imgs, varargin)
%
%  bkg = backgroundFromMin(imgs, param)
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
%
% output:
%    bkg     estimated background image 
%
% See also: backgroundFromOpening

param = parseParameter(varargin);

isTiled = false;
if isa(imgs, 'ImageSourceTiled') || isa(imgs, 'ImageSourceTagged')
   isTiled = true;
   n = obj.ncells;
elseif isa(imgs, 'ImageSource')
   imgs = imgs.celldata;
   n = numel(imgs);
else
   n = numel(imgs);
end

if n == 0
   bkg = [];
   return
end
   
if n > getParameter(param, 'maximages', 500)
   irange=randperm(n);
else
   irange = 1:n;
end

bkg = [];
for i = irange
   if isTiled
      img = imgs.data(i);
   else
      img = imgs{i};
   end

   if isempty(bkg)
      bkg = img;
   else
      bkg = min(bkg,img);
   end
end

if getParameter(param, 'filter', true)
   bkg = gaussianFilter(bkg, getParameter(param, 'filtersize', 100));
end

end
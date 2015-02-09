function img = stitchImagesByPyramid(imgs, ipos, varargin)
%
% img = stitchImagesByPyramid(imgs, shifts)
%
% description:
%     stiches images in cell array imgs together using a laplacian pyramid type blending 
%
% input:
%     imgs         images as cell array
%     ipos         positions or relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%     param        parameter struct with entries
%                  .size     final image size, ipos are interpreted as shifts if no size is specified, otherwise as absolute positions in image of this size ([])
%                  .level    depth of laplacian pyramid (3)
%                  .sigma    sigma for gaussian blur of masks, filter size is 2*sigma (15)
%
% output:
%     img          stiched image
%
% See also: stitchImages, stitchImagesByMax, stitchImagesByMin, stitchImagesByOverwrite, alignImages

if numel(imgs) == 1
   img = imgs{1};
   return
end

isizes = cellfunc(@size, imgs);

param = parseParameter(varargin);

asize = getParameter(param, 'size', []);
level  = getParameter(param, 'level', 3);
blur   = getParameter(param, 'sigma', 15);

if isempty(asize)
   [ashifts, asize] = absoluteShiftsAndSize(ipos, isizes);
else
   ashifts = ipos;
end


% note: simultaneous fuse, for larger number of images implment sequential fusion
% create masks for blending
masks = repmat({zeros(asize)}, size(imgs));
imglp = cell(size(imgs));
nd = length(isizes{1});
blurh = fspecial('gauss',2*blur,blur); % feather the border
%level = 3;
ni = numel(imgs);

%background image
imgref = stitchImagesByInterpolate(imgs, ipos, varargin{:});

for i= 1:ni
   sh = ashifts{i};
   is = isizes{i};
   mids = cell(1,nd);
   for n = 1:nd
      mids{n} = (sh(n)+ 1) : (sh(n) + is(n));
   end
   mm = masks{i};
   mm(mids{:}) = 1.;
   masks{i} = mm;
 
   imgm = imgref;
   imgm(mids{:}) = imgs{i};
   
   imglp{i} = genPyr(imgm,'lap',level);
end

%normalize masks and blur
masksn = zeros(size(masks{1}));
for i = 1:ni
   %masks{i} = imfilter(masks{i},blurh,'replicate'); 
   masksn = masksn + masks{i};
end
for i = 1:ni
   masks{i} = masks{i} ./ masksn;
   masks{i} = imfilter(masks{i},blurh,'replicate'); 
end

% figure(6);
% implottiling(masks);


imgbp = cell(1,level); % the blended pyramid
for p = 1:level
   [Mp, Np] = size(imglp{1}{p});
   imgbp{p} = zeros(Mp, Np);
   for i = 1:ni
      masks{i} = imresize(masks{i},[Mp Np]);
      imgbp{p} = imgbp{p} + imglp{i}{p}.* masks{i};
   end
end
img = pyrReconstruct(imgbp);

end
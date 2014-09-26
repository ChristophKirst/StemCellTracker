function [roi, pks] = findROIsByPeakVolume(imgs, varargin)
%
% roi = findROIsByPeakVolume(img, param)
% roi = findROIsByPeakVolume(imgs, shifts, param)
% roi = findROIsByPeakVolume(isalgn, param)
%
% descritption:
%     finds shapes in an image img bz detecting peaks and determining theri alpha volumes.
%
% input:
%    img    image as numeric array
%    imgs   images as cell arrays
%    isalgn ImageSourceAligned class
%    param  paramteer struct
%           .radius    probe radius passed to alphavol (100)
%           
%           .plot      plot the result
%           other parametre as in stitchPeaks, findPeaksByHmax
%           
% output:
%    roi    ROI objects representing the ROI
%    pks    (optional) the detected and stitched peaks
%
% See also: stitchPeaks, findPeaksByHmax

%prepare
if isnumeric(imgs)
   shifts = {zeros(1, ndims(imgs))};
   imgs = {imgs};
elseif iscell(imgs)
   if nargin < 2 || ~iscell(varargin{1})
      error('findROIsByPeakVolume: expects image shifts as second argument!')
   end
   
   shifts = varargin{1};
   varargin = varargin(2:end);
end

param = parseParameter(varargin);

% find peaks
if isa(imgs, 'ImageSourceAligned')
   n = imgs.nnodes;
   pks = cell(1, n);
   for i = 1:n
      pks{i} = findPeaksByHmax(imgs.image(i), param);
   end
   
   pks = stitchPoints(pks, imgs.absoluteShiftsAndSize, imgs.imageSizes);
   
else
   n = numel(imgs);
   pks = cell(1, n);
   for i = 1:n
      pks{i} = findPeaksByHmax(imgs{i}, param);
   end
   
   pks = stitchPoints(pks, shifts, cellfunc(@size, imgs), param);
end

% get the polygons from the points
roi = points2shapes(pks,param);
nr = length(roi);


if getParameter(param, 'plot', false)
   cc=colorcube(nr);
   
   hold on;
   for ii=1:nr
      rr = roi{ii};
      %plot(rr(1,:)', rr(2,:)', '.', 'Color', cc(ii,:));
      plot(rr(1,:)', rr(2,:)', 'LineWidth',1, 'Color', cc(ii,:));
   end
end

%convert polygons to ROIs 
for i = 1:length(roi)
   roi{i} = ROIPolygon(roi{i});
end

end





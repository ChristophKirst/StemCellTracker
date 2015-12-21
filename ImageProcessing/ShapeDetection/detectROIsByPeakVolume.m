function [roi, pks] = detectROIsByPeakVolume(imgs, varargin)
%
% roi = detectROIsByPeakVolume(img, param)
% roi = detectROIsByPeakVolume(imgs, shifts, param)
% roi = detectROIsByPeakVolume(isalgn, param)
%
% descritption:
%     finds shapes in an image img by detecting peaks and determining their alpha volumes.
%
% input:
%    img    image as numeric array
%    imgs   images as cell arrays
%    is     ImageSource or Alignment class
%    param  paramteer struct
%           .radius         probe radius passed to detectAlphaVolume (100)
%           .plot           plot the result
%           .origin.images  origin of the image set
%           .origin.plot    origin in the plot
%           other parametre as in stitchPeaks, detectPeaksByHmax
%           
% output:
%    roi    ROI objects representing the ROI
%    pks    (optional) the detected and stitched peaks
%
% See also: stitchPeaks, detectPeaksByHmax

%prepare
if isnumeric(imgs)
   shifts = {zeros(1, ndims(imgs))};
   imgs = {imgs};
   
elseif iscell(imgs)
   if nargin < 2 || ~iscell(varargin{1})
      error('detectROIsByPeakVolume: expects image shifts as second argument!')
   end
   
   shifts = varargin{1};
   varargin = varargin(2:end);
end

param = parseParameter(varargin);
origin = getParameter(param, 'origin.images', []);

% find peaks
if isa(imgs, 'Alignment')
   n = imgs.nNodes;
   pks = cell(1, n);
   parfor i = 1:n
      pks{i} = detectPeaksByHmax(imgs.nodeData(i), param); %#ok<PFBNS>
   end
   
   pks = stitchPoints(pks, imgs.imageShifts, imgs.imageSizes);
   
   if isempty(origin)
      origin = imgs.origin;
   end
   
else
   n = numel(imgs);
   pks = cell(1, n);
   parfor i = 1:n
      pks{i} = detectPeaksByHmax(imgs{i}, param);
   end
   
   pks = stitchPoints(pks, shifts, cellfunc(@size, imgs), param);
   
   if isempty(origin)
      origin = [1,1];
   end
end


if size(pks,2) < 3
   roi = [];
   return
end
   
% get the polygons from the points
roi = pointsToPolygons(pks, param);
nr = length(roi);


if getParameter(param, 'plot', false)
   cc = colorcube(nr);
   
   plorigin = getParameter(param, 'origin.plot', origin);
   
   hold on;
   for ii=1:nr
      rr = polygonShift(roi{ii}, plorigin - 1);      
      %plot(rr(1,:)' + plorigin(1) - 1, rr(2,:)' + plorigin(2) - 1, 'LineWidth', 1, 'Color', cc(ii,:));
      polygonPlot(rr, 'FaceColor' , 'none', 'EdgeColor', cc(ii,:));
      plot(pks(1,:) + plorigin(1) - 1, pks(2,:) + plorigin(2) - 1, '.', 'Color', 'k');
   end
end

%convert polygons to ROIs 
for i = 1:length(roi)
   roi{i} = ROIPolygon(roi{i});
   %roi{i} = roi{i}.shift(origin -1);
end

end





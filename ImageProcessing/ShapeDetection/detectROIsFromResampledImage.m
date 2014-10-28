function imglab = detectROIsFromResampledImage(img, scale, varargin)
%
% imglab = detectROIsFromResampledImage(img, param)
%
% descritption:
%     finds shapes in an resmapled image using morphological closing and thresholding
%     rescales rois to full scale assuming the resampling scale
%
% input:
%    img    image
%    scale  resmapling scale
%    param  parameter struct
%           .threshold    threshold ([] = thresholdFirstMin(img))
%           .resize       resize factor before detection of objects ([]=false)
%           .filtersize   size of Gaussian filter ([] = no filter) 
%           .strel        structure element for morphological closing ([] = strel('disk', 20))
%           .output       'ROI' or 'label' ('label')
%           .plot         plot result
%           options for points2shape
%           
% output:
%    imglab labeled image or array of ROIPolygon objects

param = parseParameter(varargin);

imglab = img;

fs = getParameter(param, 'filtersize', []);
if ~isempty(fs)
   imglab = filterGaussian(imglab, fs);
end

se= getParameter(param, 'strel', strel('disk', 20));
if isnumeric(se)
   se = strel('disk', se);
end
imglab = imclose(imglab, se);


th = getParameter(param, 'threshold', []);
if isempty(th)
   th = thresholdFirstMin(mat2gray(imglab));
end

imglab = imglab > th;

if isequal(getParameter(param, 'output', 'ROI'), 'ROI')

%   ids = regionprops(imglab, 'PixelIdxList');
%   ids = imind2sub(size(imglab), ids);

   if ~ismatrix(imglab)
      error('detectROIsFromResampledImage: image dimension assumed to be 2!');
   end
   
   [idx, idy] = find(imglab);

   
   if ~isempty(idx)
      
      pks = [idx'; idy'];
      
      % get the polygons from the points
      shapes = points2shapes(pks, param);
      nshapes = length(shapes);
      
      if getParameter(param, 'plot', false)
         implot(img);
         cc = colorcube(nshapes);
         for i = 1:nshapes
            shp = shapes{i};
            plot(shp(1,:)', shp(2,:)', 'LineWidth', 1, 'Color', cc(ii,:))
         end
      end

      % rescale
      for i = 1:nshapes
         shapes{i} = round(shapes{i} * (1/scale));
      end
         
      %convert polygons to ROIs 
      roi(nshapes) = ROIPolygon;
      
      for i = 1:nshapes
         roi(i) = ROIPolygon(shapes{i});
      end
   end
   
   imglab = roi;
  
else

   imglab = bwlabeln(imglab);
   
   if getParameter(param, 'plot', false)
      implot(imoverlaylabel(img, imglab));
   end
   
   imglab = imresize(imglab, 1/scale);
end
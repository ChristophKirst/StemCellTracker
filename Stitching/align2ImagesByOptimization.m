function [shift, img] = align2ImagesByOptimization(img1, img2, varargin)
%
% img = align2ImagesByOptimization(img1, img2, param)
%
% description:
%     aligns images img1, img2 using an optimization method
%
% input:
%     img1, img1   images to be alligned, assumed to have same orientation
%     param        (optional)   struct with entries
%                  .metric      metric as used by imregister
%                  .optimizer   optimizer as used by imregister
%                  .overlap.max maximal overlap of images in selected direction in pixel (150)
%                  .direction   'lr', 'rl', 'bt', 'tb' 'du', 'ud' ('lr' = img1 left, img2 right)
%                               [1 0 (0)], [-1 0 (0)], [0 1 (0)], [0 -1 (0)], [0 0 1], [0 0 -1]
% 
% output:
%     shift        the shift between origin of img1 to origin of img2 in pixel coordinates and pixel units
%     img          (optional) combined image using stitch2ImagesByOverwrite
%
% See also:  align2ImagesByShift, stitch2ImagesByMean, stich2ImagesByWatershed, alignImages, stitchImages, imregister

param = parseParameter(varargin{:});

si1 = size(img1);
si2 = size(img2);

dim = length(si1);
if dim ~= length(si2)
   error('alignImages: images dimension mistmatch!');
end

maxovl = getParameter(param, {'overlap', 'max'}, 150);
if isempty(maxovl) 
   maxovl = 150;
end
maxovl = maxovl(1);
%maxovl = maxovl(:)';
%maxovl = padarray(maxovl, [0, max(0, dim - length(maxovl))], 'circular', 'post');

% minovl = getParameter(param, {'overlap', 'min'}, 10);
% if isempty(minovl) 
%    minovl = 10;
% end
% minovl = minovl(:)';
% minovl = padarray(minovl, [0, max(0, dim - length(minovl))], 'circular', 'post');

direct = getParameter(param, 'direction', 'lr');
if isempty(direct)
   direct = 'lr';
end

if ischar(direct)
   switch dim
      case 2
         if  ~any(strcmp({'lr', 'rl', 'tb', 'bt'}, direct))
            error('align2Images: direction %s not valid');
         end
      case 3
         if  ~any(strcmp({'lr', 'rl', 'tb', 'bt', 'ud', 'du'}, direct))
            error('align2Images: direction %s not valid');
         end
   end
   
   switch direct
      case 'lr'
         direct = 1;
      case 'rl'
         direct = -1;
      case 'bt'
         direct = 2;
      case 'tb'
         direct = -2;
      case 'du'
         direct = 3;
      case 'ud'
         direct = -3;    
   end
elseif isnumeric(direct) && length(direct) > 1
   p = find(direct);
   if length(p) > 1
      error('align2Images: direction not valid');
   end
   direct = p * sign(direct(p));
end

if ~isscalar(direct) ||  direct < 1 || direct > dim
   error('align2Images: direction not valid'); 
end


%[optimizer, metric] = imregconfig('multimodal');
[optimizer, metric] = imregconfig('monomodal');

%%% for monomodal
optimizer.MaximumIterations = 3000;
optimizer.RelaxationFactor = 0.95;
optimizer.MaximumStepLength = 0.5;
optimizer.MinimumStepLength = 0.05;

%%% for gradient descent
%optimizer.InitialRadius = 0.009;
%optimizer.Epsilon = 1.5e-4;
%optimizer.GrowthFactor = 1.01;
%optimizer.GradientMagnitudeTolerance = 100;

%%% for evolutionary algorithm
%optimizer.InitialRadius = 15;
%optimizer.MaximumIterations = 3000;
%metric.NumberOfHistogramBins = 256;
%metric.NumberOfSpatialSamples = 1000;

poptimizer = getParameter(param, 'optimizer', optimizer);
pmetric = getParameter(param, 'metric', metric);

if ~isempty(poptimizer)
   optimizer = poptimizer;
end
if ~isempty(pmetric)
   metric = pmetric;
end






% extract the borders

switch dim
   case 2
      switch direct
         case 1
            img1b = img1(end-maxovl(1)+1:end, :);
            img2b = img2(1:maxovl(1), :);
         case -1
            img2b = img2(end-maxovl(1)+1:end,:);
            img1b = img1(1:maxovl(1),:);
         case 2
            img1b = img1(:, end-maxovl(1)+1:end);
            img2b = img2(:, 1:maxovl(1));
         case -2
            img2b = img2(:, end-maxovl(1)+1:end);
            img1b = img1(:, 1:maxovl(1));
         otherwise
            error('align2Images: direction not valid!');
      end
      
   case 3
      switch direct
         case 1
            img1b = img1(end-maxovl(1)+1:end, :, :);
            img2b = img2(1:maxovl(1), :, :);
         case -1
            img2b = img2(end-maxovl(1)+1:end,:, :);
            img1b = img1(1:maxovl(1),:, :);
         case 2
            img1b = img1(:, end-maxovl(1)+1:end, :);
            img2b = img2(:, 1:maxovl(1), :);
         case -2
            img2b = img2(:, end-maxovl(1)+1:end, :);
            img1b = img1(:, 1:maxovl(1), :);
         case 3
            img1b = img1(:, :, end-maxovl(1)+1:end, :);
            img2b = img2(:, :, 1:maxovl(1));
         case -3
            img2b = img2(:, end-maxovl(1)+1:end, :);
            img1b = img1(:, :, 1:maxovl(1));
         otherwise
            error('align2Images: direction not valid!');
      end
end   
   
% figure(13)
% clf
% imsubplot(1,2,1)
% implot(img1b)
% imsubplot(1,2,2);
% implot(img2b)

% optimize - note: the if max overlap is choosen to big the alorith might not converge !
%tform = imregtform(img2b, img1b, 'translation', optimizer, metric, 'InitialTransform', trsf0,'DisplayOptimization',false)
tform = imregtform(img2b, img1b, 'translation', optimizer, metric,'DisplayOptimization',false);


%calculate shifts for full images
shift = tform.T;

switch dim
   case 2
      shift = shift(3,[2 1]);
   case 3
      shift = shift(4,[2 1 3]);
end

%pixel shifts
shift = round(shift);

%absolute shift due to orientation and imgage size
direct = abs(direct);
ev = zeros(1,dim); ev(direct) = 1;

if sum(direct) > 0
   shift = shift + (size(img1,direct) - (maxovl )) * ev;
else
   shift = shift - (size(img2,direct) + (maxovl)) * ev;
end


if nargout > 1
   img = stitch2ImagesByOverwrite(img1, img2, shift);
end
   
end
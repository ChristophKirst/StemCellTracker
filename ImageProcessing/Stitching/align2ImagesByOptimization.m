function [shift, quality] = align2ImagesByOptimization(img1, img2, varargin)
%
% [shift, quality] = align2ImagesByOptimization(img1, img2, param)
%
% description:
%     finds best global alignment of two full images img1 and img2 
%
% input:
%     img1, img1   images to be alligned, img1 is left img 2 is right
%     param        (optional)   struct with entries
%                  .metric      metric as used by imregister
%                  .optimizer   optimizer as used by imregister
%
% output:
%     shift        the shift between origin of img1 to origin of img2 in pixel coordinates and pixel units
%     quality      (optional) quality measure of the alignment = 1
%
% See also:  align2ImagesByOptimization, imregister

param = parseParameter(varargin{:});

si1 = size(img1);
si2 = size(img2);

dim = length(si1);
if dim ~= length(si2)
   error('alignByOptimization: image dimension mistmatch, %g ~= %g!', dim, length(si2));
end

%[optimizer, metric] = imregconfig('multimodal');
[optimizer, metric] = imregconfig('monomodal'); % images are from the same microscope 

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

%tform = imregtform(img2b, img1b, 'translation', optimizer, metric, 'InitialTransform', trsf0,'DisplayOptimization',false)
tform = imregtform(img2, img1, 'translation', optimizer, metric,'DisplayOptimization',false);

%calculate shifts for full images
shift = tform.T;

switch dim
   case 2
      shift = shift(3,[2 1]);
   case 3
      shift = shift(4,[2 1 3]);
end

shift = round(shift);
quality = 1;

end
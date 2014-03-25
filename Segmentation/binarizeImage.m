function bw = binarizeImage(image, threshold, varargin)
%
% bw = binarizeImage(image, threshold)
%
% description: 
%     binarize image using method or level threshold
%
% input:
%    image      image to be binarized
%    threshold  method string or numeric threshold level
%               method can be 'Otsu', 'Entropy', 'MixtureOfGaussians', 'MutualEntropy'
%
% output:
%     bw        binarized image
%
% See also: im2bw, otsuThreshold, entropyThreshold, mixtureOfGaussiansThreshold, mutualEntropyThreshold

if nargin < 1
   threshold = 'Otsu';
end

if ~isnumeric(threshold)
   if sum(strcmp({'Otsu' 'Entropy' 'MixtureOfGaussians', 'MutualEntropy'}, threshold)) > 0
      methodName = [lower(threshold(1)) threshold(2:end) 'Threshold'];
      thresFunc = str2func(methodName);
      threshold = thresFunc(image, varargin{:});      
   else
      error(['binarizeImage: no available method ' threshold ' for thresholding!'])
   end
end

bw = im2bw(image, threshold);

end

   
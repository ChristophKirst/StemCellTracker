function threshold = mixtureOfGaussiansThreshold(image, probability)
%
% threshold = mixtureOfGaussiansThreshold(image, probability)
%
% description:
%      fits the sum of three gaussians to the histogram of image and
%      infers the threshold based on the foreground probablity
%
% input: 
%      image        input image to calculate threshold for
%      probablity   probability of a pixel belonging to foreground
%
% output:
%      threshold    threshold level
%
% reference:
%      CellProfiler, Broad Institute, MIT

if nargin < 2
   probability = 0.5;
end
probBackground = 1 - probability;

if ~ismatrix(image)
   image = rgb2gray(image);
end
image = double(image(:));

if length(image) > 512^2
   %defaultStream = RandStream.getGlobalStream;
   %savedState = defaultStream.State;
   %RandStream.setGlobalStream(RandStream('mt19937ar','seed',0));
   indexes = randperm(length(image));
   %defaultStream.State = savedState;
   image = image(indexes(1:512^2));
end

image = sort(image);
image = image(ceil(length(image)*0.01):round(length(image)*0.99)); % cut off lower and top 1% 


nClasses = 3;
classMean(1) = image(round(length(image)*probBackground/2));                      % Initialize background class
classMean(3) = image(round(length(image)*(1 - probability/2)));                % Initialize foregournd class
classMean(2) = (classMean(1) + classMean(3))/2;                                   % Initialize intermediate class
classStd(1:3) = 0.15;

pClass(1) = 3/4*probBackground;
pClass(2) = 1/4*probBackground + 1/4*probability;
pClass(3) = 3/4*probability;


 % Expectation-Maximization algorithm for fitting the three Gaussians
 delta = 1;
 while delta > 0.001
     oldClassMean = classMean;

     % Update probabilities of a pixel belonging to the background or
     % object1 or object2
     for k = 1:nClasses
         pPixelClass(:,k) = pClass(k)* 1/sqrt(2*pi*classStd(k)^2) * exp(-(image - classMean(k)).^2/(2*classStd(k)^2)); %#ok<AGROW>
     end
     pPixelClass = pPixelClass ./ repmat(sum(pPixelClass,2) + eps,[1 nClasses]);

     % Update parameters in Gaussian distributions
     for k = 1:nClasses
         pClass(k) = mean(pPixelClass(:,k)); %#ok<AGROW>
         classMean(k) = sum(pPixelClass(:,k).*image)/(length(image)*pClass(k)); %#ok<AGROW>
         classStd(k)  = sqrt(sum(pPixelClass(:,k).*(image - classMean(k)).^2)/(length(image)*pClass(k))) + sqrt(eps);    % Add sqrt(eps) to avoid division by zero
     end

     % Calculate change
     delta = sum(abs(classMean - oldClassMean));
 end

 % find if intermediate class 2 encodes background or object pixels.
 threshold = linspace(classMean(1),classMean(3),10000);
 class1Gaussian = pClass(1) * 1/sqrt(2*pi*classStd(1)^2) * exp(-(threshold - classMean(1)).^2/(2*classStd(1)^2));
 class2Gaussian = pClass(2) * 1/sqrt(2*pi*classStd(2)^2) * exp(-(threshold - classMean(2)).^2/(2*classStd(2)^2));
 class3Gaussian = pClass(3) * 1/sqrt(2*pi*classStd(3)^2) * exp(-(threshold - classMean(3)).^2/(2*classStd(3)^2));
 if abs(pClass(2) + pClass(3) - probability) < abs(pClass(3) - probability)
     % intermediate class 2 encodes object pixels
     backgroundDistribution = class1Gaussian;
     objectDistribution = class2Gaussian + class3Gaussian;
 else
     % Intermediate class 2 encodes background pixels
     backgroundDistribution = class1Gaussian + class2Gaussian;
     objectDistribution = class3Gaussian;
 end

 % threshold at the intersection of the background and foregournd
 [ignore,index] = min(abs(backgroundDistribution - objectDistribution)); %#ok Ignore MLint
 threshold = threshold(index);
 
end


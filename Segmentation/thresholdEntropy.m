function threshold = thresholdEntropy(image, param)
%
% threshold = thresholdEntropy(image, param)
%
% description:
%    finds a threshold based on entropies of the image histrogram
%
% input:
%    image       the image to calculate threshold for
%    param       paramter struct with entry
%                .minlog2val (-10)  or the value directly
%
% output:
%    threshold   calculated threshold
%
% reference:
%    Kapur, Sahoo, & Wong, A new Method for Gray Level Picture Thresholding
%    Using the Entropy of the Histogram, 1985

if nargin < 2 
   param = [];
end

if ~isnumeric(param) || isempty(param)
   minlog2 = getParameter(param, 'minlog2val', -10);
else
   minlog2 = param;
end

image = double(image);

%histrogram
log2img = log2(image(:)+eps);
log2img(log2img < minlog2) = minlog2;
[N, X] = hist(log2img, 256);

% drop any zero bins
drop = (N == 0);
N(drop) = [];
X(drop) = [];

% check for corner cases
if length(X) == 1
    threshold = X(1);
    return;
end

% Normalize to probabilities
P = N / sum(N);

% Find the probabilities totals up to and above each possible threshold.
loSum = cumsum(P);
hiSum = loSum(end) - loSum;
loE = cumsum(P .* log2(P));
hiE = loE(end) - loE;

% compute the entropies
loEntropy = loE ./ loSum - log2(loSum);
hiEntropy = hiE ./ hiSum - log2(hiSum);

sumEntropy = loEntropy(1:end-1) + hiEntropy(1:end-1);
sumEntropy(~ isfinite(sumEntropy)) = Inf;
entry = find(sumEntropy == min(sumEntropy), 1, 'first');
threshold = 2^((X(entry) + X(entry+1)) / 2);

end










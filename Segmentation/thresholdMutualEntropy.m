function threshold = thresholdMutualEntropy(image)
%
% threshold = thresholdMutualEntropy( image )
%
% description:
%    finds a threshold based on entropies of the image histrogram
%
% input:
%    image       the image to calculate threshold for
%
% output:
%    threshold   calculated threshold
%
% reference:
%    Johannsen, Bille, A threshold selection method using information measures
%    Proceedings 6th Int. conf. Pattern Recognition, 1982

if nargin < 2 
   param = [];
end

if ~isnumeric(param)
   minlog2 = getParameter(param, 'minlog2val', -10);
else
   minlog2 = param;
end
   

image = double(image);

%histrogram
log2img = log2(image(:)+eps);
log2img(log2img < minlog2) = minlog2;
[N, X] = hist(log2img, 256);


%histrogram
[N, X] = hist(log2(image(:)), 256);

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
loPSum = cumsum(P);
hiPSum = loPSum(end) - loPSum;


%calculate loE and hiE for each threshold t

loS = log(loPSum) - (P .* log(P) + loPSum .* log(loPSum)) ./ loPSum;
hiS = log(hiPSum) - (P .* log(P) + hiPSum .* log(hiPSum)) ./ hiPSum;


sumS = loS + hiS;
sumS(~ isfinite(sumS)) = Inf;

entry = find(sumS == min(sumS), 1, 'first');
threshold = 2^X(entry);

end


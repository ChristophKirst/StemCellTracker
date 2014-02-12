function threshold = mutualEntropyThreshold(image)
%
% threshold = mutualEntropyThreshold( image )
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


if ~ismatrix(image)
   image = rgb2gray(image);
end
image = double(image);

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


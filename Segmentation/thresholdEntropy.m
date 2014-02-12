function threshold = entropyThreshold(image)
%
% threshold = entropyThreshold( image )
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
%    Kapur, Sahoo, & Wong, A new Method for Gray Level Picture Thresholding
%    Using the Entropy of the Histogram, 1985


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










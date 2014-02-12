function y = randsample(n, k)
%
% y = randsample(n, k)
%
% descritption:
%    draws k indices without repititions from 1:n
%
% input:
%    n    number of elements
%    k    number of elemnts to draw
%
% output:
%    y    random sample


% If the sample is a sizeable fraction of the population, just
% randomize the whole population (which involves a full sort
% of n random values), and take the first k.
if 4*k > n
    rp = randperm(n);
    y = rp(1:k);

% If the sample is a small fraction of the population, a full
% sort is wasteful.  Repeatedly sample with replacement until
% there are k unique values.
else
    x = zeros(1,n); % flags
    sumx = 0;
    while sumx < k
        x(ceil(n * rand(1,k-sumx))) = 1; % sample w/replacement
        sumx = sum(x); % count how many unique elements so far
    end
    y = find(x > 0);
    y = y(randperm(k));
end

% a scalar loop version
%
% x = 1:n;
% n = n:(-1):(n-k+1);
% y = zeros(1,k);
% j = ceil(n .* rand(1,k));
% for i = 1:k
%     y(i) = x(j(i));
%     x(j(i)) = x(n(i));
% end

end
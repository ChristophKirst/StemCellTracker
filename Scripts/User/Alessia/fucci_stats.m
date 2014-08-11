function stats = fucci_stats(imglab, rcenter, stats, varargin)

if nargin > 1
   verbose = varargin{1};
else
   verbose = false;
end
   

%% Distance from Center

%rcenter = [1100, 555];

stats = imstatistics(imglab, stats, 'Centroid');
pos = [stats.Centroid];

for l=length(stats):-1:1 
   vec = pos(:,l) - rcenter';
   dist(l) = norm(vec);
   theta(l) = atan2(vec(2), vec(1));
   if theta(l) < 0
      theta(l) = theta(l) + 2 * pi;
   end
end

dist = num2cell(dist);
[stats.dist] = dist{:};

theta = num2cell(theta);
[stats.theta] = theta{:};


end


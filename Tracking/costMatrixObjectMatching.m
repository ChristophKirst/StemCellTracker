function cost = costMatrixObjectMatching(data0, data1, param)
% 
% cost = costMatrixObjectMatching(data0, data1, param)
%
% description:
%   calculates the cost for matching the objects in data0 to data1
%
% input:
%   data*        classes that contain entry objects or directly
%                arrays of objects to calculate linking cost 
%                objects used entries: r, time, volume, intensity, identity
%   param        parameter structure with entries:
%                .cost.creation           - costs for creating new objects ([])
%                .cost.deletion           - costs for deleting objects ([]) 
%                .cost.join               - cost for join objects (Inf = not allowed)
%                .cost.split              - cost for splitting objects (Inf = not allowed)
%                .cutoff.cost             - overall cost cutoff ([] = automatic)
%                .cutoff.dist             - cutoff for spatial distance ([] = automatic)
%                .cutoff.time             - cutoff for temporal distances (Inf)
%                .cutoff.join.distance    - cutoff for spatial distance when joining objects ([] = automatic)
%                .cutoff.split.distance   - cutoff for spatial distance when joining objects ([] = automatic)
%                .cutoff.join.nmax        - cutoff for maximal number of objects to join (2)
%                .cutoff.split.nmax       - cutoff for maximal number of objects to split into (2)
%                .weight.dist             - weight for distances in space (1.0)
%                .weight.time.forward     - weight for positive time distances (0.0)
%                .weight.time.backward    - weight for negative time distances (0.0)
%                .weight.volume           - weight for distances in volume (0.0)
%                .weight.intensity        - weight for distances in intensity (0.0)
%                .weight.type             - weight for different types of objects (Inf)
%                .weight.join.volume      - weight for volume differnces of joining objects
%                .weight.join.intensity   - weight for intensity differences for joining objects
%                .weight.join.distance    - weight for distance differences for joining objects
%                .weight.split.volume     - weight for volume differnces of splitting objects
%                .weight.split.intensity  - weight for intensity differences for splitting objects
%                .weight.split.distance   - weight for distance differences for splitting objects
%
% output:
%   cost         the cost matrix
%
% todo:
%   motion propagation / creation and deletion cost using distance to the image borders
%   optimize by reducing calcualtions to groups of same identity
%
% See also: distanceMatrix, matchObjects, estimateDistanceCutoff, estimateNonMatchingCost

% initialize

if nargin < 3
   param = [];
end

% basic join of objects
dist_cutoff = getParameter(param, {'cutoff', 'dist'}, []);      % cutoff for weighted distances in space
time_cutoff = getParameter(param, {'cutoff', 'time'}, Inf);     % cutoff for weighted distances in time
cost_cutoff = getParameter(param, {'cutoff', 'cost'}, []);      % over all cost cutoff 

cost_creation = getParameter(param, {'cost', 'creation'}, []);  % cost of creating objects
cost_deletion = getParameter(param, {'cost', 'deletion'}, []);  % cost of deleting objects

w_dist          = getParameter(param, {'weight', 'dist'}, 1.0);              % weight for spatial distance
w_time_forward  = getParameter(param, {'weight', 'time', 'forward'}, Inf);   % weight for forward time differences
w_time_backward = getParameter(param, {'weight', 'time', 'backward'}, Inf);  % weight for backward time differences
w_volume        = getParameter(param, {'weight', 'volume'}, 0.0);            % weight for size distance
w_intensity     = getParameter(param, {'weight', 'intensity'}, 0.0);         % weight for intensity distance
w_type          = getParameter(param, {'weight', 'type'}, Inf);              % weight for two objects that have different identity


%joining an splitting
join_cost = getParameter(param, {'cost', 'join'}, Inf); 
split_cost = getParameter(param, {'cost', 'split'}, Inf);

join_nmax = getParameter(param, {'cutoff', 'join', 'namx'}, 2);
split_namx = getParameter(param, {'cutoff', 'split', 'namx'}, 2);

join_cutoff_distance = getParameter(param, {'cutoff', 'join', 'distance'}, []);
split_cutoff_distance = getParameter(param, {'cutoff', 'split', 'distnace'}, []);

join_w_distance = getParameter(param, {'weight', 'join', 'distance'}, 1.0);
split_w_distance = getParameter(param, {'weight', 'split', 'distance'}, 1.0);
join_w_volume = getParameter(param, {'weight', 'join', 'volume'}, 0.0);
split_w_volume = getParameter(param, {'weight', 'split', 'volume'}, 0.0);
join_w_intensity = getParameter(param, {'weight', 'join', 'intensity'}, 0.0);
split_w_intensity = getParameter(param, {'weight', 'split', 'intensity'}, 0.0);



if isentry(data0, 'objects')
   data0 = data0.objects;
end
if isentry(data1, 'objects')
   data1 = data1.objects;
end

n0 = length(data0);
n1 = length(data1);


% calculate spatial distances 

if isentry(data0, 'r') && ( w_dist > 0 )
   r0 = [ data0.r ];
   r1 = [ data1.r ];
   
   if size(r0,2) ~= n0 || size(r1,2) ~= n1
      error('positions in data do not have correct dimensions: %d ~= %d or %d =~ %d', size(r0,2), n0, size(r1,2), n1);
   end
  
   dim = data0.dim;
   if dim ~= data1.dim
      error('positions in data have different dimensions: %d ~= %d', dim, data1.dim);
   end
   
   cost = w_dist * sqrt(sum(bsxfun(@minus,reshape(r0',[n0,1,dim]),reshape(r1',[1,n1,dim])).^2,3));

   if isempty(dist_cutoff)
      dist_cutoff = estimateDistanceCutoff(cost, param);
   end
   if dist_cutoff < Inf
      cost(cost > dist_cutoff) = Inf; 
   end   
else
   cost = zeros(n0, n1);
end


% calculate cost for time travel

if isfield(data0, 'time') && (w_time_forward > 0 || w_time_backward > 0)
   
   indx = find(cost < Inf);
   [row, col] = ind2sub(size(cost), indx);
  
   t0 = [ data0.time ];
   t1 = [ data1.time ];
   
   if length(t0) ~= n0 || length(t1) ~= n1
      error('times in data do not have correct dimensions: %d ~= %d or %d =~ %d', length(t0), n0, length(t1), n1);
   end
   
   dt = t0(row) - t1(col); 
   ind = fnd(dt > 0 );
   dt(ind)  = w_time_forward * dt(ind);
   
   ind = find(dt < 0 );
   dt(ind) = w_time_backward * abs(dt(ind));
   
   if ~isempty(time_cutoff) 
      dt(dt > time_cutoff) = Inf;
   end
   
   cost(indx) = cost(indx) + dt;
end


% calculate cost for other differences

if (w_volume > 0 || w_intensity > 0 || w_type > 0)
   indx = find(cost < Inf);
   [row, col] = ind2sub(size(cost), indx);
end

if isfield(data0, 'volume') && w_volume > 0 
   s0 = [ data0.size ];
   s1 = [ data1.size ];

   if length(s0) ~= n0 || length(s1) ~= n1
      error('sizes in data do not have correct dimensions: %d ~= %d or %d =~ %d', length(s0), n0, length(s1), n1);
   end

   cost(indx) = cost(indx) + w_volume * abs( s0(row) - s1(col) );     
end

if isfield(data0, 'intensity') && w_intensity > 0    
   i0 = [ data0.intenisty ];
   i1 = [ data1.intensity ];

   if length(i0) ~= n0 || length(i1) ~= n1
      error('intensities in data do not have correct dimensions: %d ~= %d or %d =~ %d', length(i0), n0, length(i1), n1);
   end

   cost(indx) = cost(indx) + w_intensity * abs( i0(row) - i1(col) );    
end

if isfield(data0, 'type') && w_type > 0  
      t0 = [ data0.type ];
      t1 = [ data1.type ];
      
      if size(t0,2) ~= n0 || size(t1,2) ~= n1
         error('identities in data do not have correct dimensions: %d ~= %d or %d =~ %d', size(i0,2), n0, size(i1,2), n1);
      end

      if size(t0,1) ~= size(t1,1)
         error('identities in data have different dimensions: %d ~= %d', size(i0,1), size(i1,1));
      end

      di = sum( abs(i0(row,:) - i1(col, :)) );
      di(di > 0) = w_type;
      cost(indx) = cost(indx) + di;
end


%threshold

if ~isempty(cost_cutoff)
   cost(cost > cost_cutoff) = Inf;
end


% object creation and deletion

if isempty(cost_creation)
   cost_creation = estimateNonMatchingCost(cost, param);
end
if isempty(cost_deletion)
   cost_creation = cost_creation;
end

cost(:,n1+1) = cost_creation * ones(n0,1);

cost(n0+1,:) = cost_creation * ones(1,n1+1);

cost(n0+1,n1+1) = Inf;


% joning and splitting cost






end

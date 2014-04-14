function imgseg = segmentOnTiling(image, segmenter, param)
%
% imgseg = segmentOnTiling(image, segmenter, param)
%
% description:
%    uses the segmentation routine segmenter to segment on a finer tiling 
%    for paralleization and or smaller memmories
%    
% input:
%    image     3d grayscale or color image
%    segmenter function applied to image and param that retunrs a labeled image
%    param     (optional) parameter struct with entries
%              .tiles           [np, nq (, nl)]  number of tiles in the dimension
%              .objects.size     [dp, dq, dl] expected size of objects to estimate sensible values for all parameters below
%              .overlap         overlapp in pixel [op,pq(,ol)]
%              .join.offset     centroid offset from border for objects to be considered for joining
%              .join.distance   maximal spatial distance for two centroids of objects to be considred for joining
%              .join.error      maximal number of non-overlapping pixels for objects to be join
%              .join.assing     how to assing undecided overlapp: 'random', 'min' (to smaller object) , 'max' (to bigger object)
%   
% output:
%    imgseg    full segmented image


% parameter
if nargin < 3
   param = [];
end

isize = size(image);
idim = length(isize);

obj_size = getParameter(param, {'objects', 'size'}, []);


% guess parameter based on object size
if ~isempty(obj_size)
  
   obj_size = obj_size(:);
   obj_size = [obj_size, isize((length(obj_size)+1):idim)];
   
   overlap = 5 * obj_size;
   
   join_offset = 1.5 * obj_size;
   join_distance = 1.5 * obj_size;
   join_error = 0.25 * prod(obj_size)/2 * 4 /3 * pi; % 25% ofthe volume of the ellipsoid of objsize

else   
   overlap = min(max(10, isize/10), 100);
   join_offset = round(overlap / 3);
   join_distance = round(overlap / 3);
   join_error = 0.25 * prod(overlap)/2 * 4 /3 * pi;
end

     
tiles         = getParameter(param, {'tiles'}, 1); 
overlap       = getParameter(param, {'overlap'}, overlap);
join_offset   = getParameter(param, {'join', 'offset'},join_offset);
join_distance = getParameter(param, {'join', 'offset'},join_distance);
join_error    = getParameter(param, {'join', 'error'},join_error);
join_assing   = getParameter(param, {'join', 'assing'}, 'min');


% generate dimensions for subimages and neighbour relations
tiles = [tiles(:)', ones(1, idim - length(tiles))];
   
for d = 1:idim
   hi{d} = cumsum(subsetsizes(isize(d), tiles(d))); 
   lo{d} = [1, 1 + hi{d}(1:(tiles(d)-1))];
   hi{d} = min(hi{d} + overlap(d), isize(d));
   %lo{d} = max(lo{d} - overlap(d), 1);
end

if idim == 2

   k = 1;
   p = 1;
   for i = 1:tiles(1)
      for j = 1:tiles(2)
         ijmap(i,j) = k;        
         ssize{k} = [lo{1}(i), lo{2}(j), hi{1}(i), hi{2}(j)];
         
         if i  > 1
            pairs{p} = [ijmap(i-1,j), k, 1]; %pair id 1 paird id 2, dim
            p = p + 1;
         
         end
      
         if j > 1
            pairs{p} = [ijmap(i, j-1), k, 2];
            p = p + 1;
         end
         
         k = k+1;
      end
   end
   
elseif idim == 3
   
   k = 1;
   p = 1;
   for i = 1:tiles(1)
      for j = 1:tiles(2)
         for l = 1:tiles(3)
            ijmap(i,j,l) = k;        
            ssize{k} = [lo{1}(i), lo{2}(j), lo{3}(l), hi{1}(i), hi{2}(j), hi{3}(l)];
         
            if i  > 1
               pairs{p} = [ijmap(i-1,j, l), k, 1];
               p = p + 1;
            end
      
            if j > 1
               pairs{p} = [ijmap(i, j-1, l), k, 2];
               p = p + 1;
            end
            
                  
            if l > 1
               pairs{p} = [ijmap(i, j, l-1), k, 3];
               p = p + 1;
            end
         
            k = k+1;
         end
      end
   end

else 
   error('segmentOnTiling: image dimension not 2 or 3');  
end

ntiles = k-1;
npairs = p-1;


% run segmenter on each subimage and calculate statistics

for i = 1:ntiles
   ims{i} = segmenter(imextract(img, ssize{i}), param);
   stats{i} = imstatistics(ims{i},{'PixelIdxList', 'Centroids'});
end


% join the lables

maxlabel = 1;

for p = 1:npairs
   p1 = pairs{p}(1); p2 = pairs{p}(2);
   d1 = pairs{p}(3); d2 = pairs{p}(3);
   is1 = size(ims{p1}); is2 = size(ims{p2});

   % find overlapping regoin
   rl1 = ones(1,idim); rh1 = is1; rh1(d1) = is1(d1) - overlap(d1) + 1;      
   rl2 = ones(1,idim); rh2 = is2; rh2(d2) = overlap(d2);
   
   % find labels in overlapping region
   c1 = [stats{p1}.Centroid]; 
   c2 = [stats{p2}.Centroid];
    
   id1 = ones(1,idim); id2 = id1;
   for d = 1:idim
      id1 = id1 & (c1(d,:) >= rl1(d)) & (c1(d,:) <= rh1(d));
      id2 = id2 & (c2(d,:) >= rl2(d)) & (c2(d,:) <= rh2(d));
   end
   
   c1 = c1(:, id1); c2 = c2(:, id2);
   id1 = find(id1); id2 = find(id2);
   
   % calculate distnaces
   if join_distance
   
   
   
   
end






end

function split = subsetsizes(n, p)
   q = ceil(n/p);
   l = q*p - n;
   split = q * ones(1,p);
   split(1:l) = q -1;
end


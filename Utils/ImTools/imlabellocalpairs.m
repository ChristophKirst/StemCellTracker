function [pairs, distances] = imlabellocalpairs(label, ksize, method)
%
% pairs = imlabellocalpairs(label, ksize, method)
%
% description:
%      determines all pairs of label indices that are within a box of size ksize
%      centered at the label pixel.
%
% input:
%      label     labeled image
%      ksize     p x q ( x l) box to search for neighbours, or maximal distance when using method 'distance'
%      method    (optional) either 'box' or 'distance' ('box')
% 
% output:
%      pairs     array of the form [lab1a, lab1b; lab2a, lab2b] of all found edges
%      distances (optional) the distance between each of the pairs if 'distance' method is used
%
% See also: distanceMatrix, imextractbox

if nargin < 3
   method = 'box';
end

switch method
   case 'box'
      pairs = pairs_box(label, ksize);
   otherwise
      [pairs, distances] = pairs_dist(label, ksize);
end

% alternative: via graycomatrix

end
      


function pairs = pairs_box(label, ksize)

   lbs = imlabel(label);
   pos = immask2coords(label);

   pairs = [];
   for i = 1:size(pos,1);
      box = imextractbox(label, pos(i,:), ksize);
      idx = unique(box(box > 0));
      pairs = vercat(pairs, [repmat(lbs(i), length(idx), 1), idx(:)] );
      idx = num2cell(pos(i,:));
      label(idx{:}) = 0; % delete label as we have found all neighbours
   end  
end


function [pairs, dist] = pairs_dist(label, ksize)
   dist = imlabeldistances(label);
   [i, j] = find(dist < ksize & dist > 0);   
   pairs = sort([i, j],2);
   pairs = unique(pairs, 'rows');
   
   %lbs = imlabel(label);
   %pairs = lbs(pairs) 
   dist = dist(sub2ind(size(dist), pairs(:,1), pairs(:,2)));
   
end
  
   

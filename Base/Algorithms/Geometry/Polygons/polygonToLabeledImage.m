function imgmask = polygonToLabeledImage(pol, varargin)
%
% imgmask = polygonToLabeledImage(pol, varargin)
%
% description:
%     converts a polygon into a mask
%
% input:
%     pol     cell array of polygons
%     param   parameter struct with entries
%             size     image size ([] = bounding box of the poly)
%             range    coords of the image as [lower, upper] ([] = bounding box of the poly)
%             
% output:
%     imglab  labled image
%
% See also: polygonFromLabeledImage, polygonFromMask

param = parseParameter(varargin);

si = getParameter(param, 'size', []);
rg = getParameter(param, 'range', []);

if isempty(rg)
   bb = polygonToBoundingBox(pol);
   rrg = sort(bb(:, [1,3]),2); 
end

dr = diff(rrg');

if isempty(si)
   si = diff(rrg');
   si(si < 1) = 1;
else
   if isempty(rg)
      rg = [0.5 * [1;1], si(:) + 0.5];
      dr = diff(rg');
   else
      rg = rrg;
   end
end
si = si(:);

if isempty(pol)
   imgmask = zeros(si');
   return
end

% scale factor
% range rx -> rx2 is mapped to 0.5 -> sx + 0.5
sc = si ./ dr(:);

% simplify, get tree and scale
nlab = length(pol);
ppol = cell(1, nlab);
tree = cell(1, nlab);

for i = 1:nlab
   [ppol{i}, tree{i}] = polygonSimplify(pol{i});
   ppol{i} = cellfunc(@(x) (x - repmat(rg(:,1), 1, size(x,2))) .* repmat(sc, 1, size(x,2)) + 0.5, ppol{i});
end

pMaskFcn = @(xy)poly2mask(xy(1,:)',xy(2,:), si(1), si(2));

imgmask = zeros(si');

for i = 1:nlab

   npol = length(ppol{i});
   if ~isempty(ppol{i})
    
      treei = logical(tree{i});
      poli = ppol{i};
   
      % build image by
      [rr,~] = find(treei);
      listToAdd = false(1,npol);
      listToAdd(setdiff(1:npol,rr)) = true;

      msk = false(si');
      while any(listToAdd)
         % Add the next enclosing mask
         nextAdd = find(listToAdd,1);
         msk = pMaskFcn(poli{nextAdd});         
         listToAdd(nextAdd) = false;
         % Subtract any enclosed masks (and then add their children)
         for nextSub = find(treei(:,nextAdd))'
            msk = msk & ~pMaskFcn(poli{nextSub});
            % Add any children
            listToAdd(treei(:, nextSub)) = true;
         end
      end
      
      imgmask(msk) = i;
   end
   
end

end



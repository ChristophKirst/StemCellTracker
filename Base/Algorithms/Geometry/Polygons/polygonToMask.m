function imgmask = polygonToMask(pol, varargin)
%
% imgmask = polygonToMask(pol, varargin)
%
% description:
%     converts a polygon into a mask
%
% input:
%     pol     polygon
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

% simplify and get tree
[pol, tree] = polygonSimplify(pol);


% scale the polygons appropiately 
% range rx -> rx2 is mapped to 0.5 sx + 0.5
sc = si ./ dr(:); 
pol = cellfunc(@(x) (x - repmat(rg(:,1), 1, size(x,2))) .* repmat(sc, 1, size(x,2)) + 0.5, pol);

pMaskFcn = @(xy)poly2mask(xy(1,:)',xy(2,:), si(1), si(2));

imgmask = false(si');
npol = length(pol);

% build image by 
[rr,~] = find(tree);
listToAdd = false(1,npol);
listToAdd(setdiff(1:npol,rr)) = true;

tree = logical(tree);

while any(listToAdd)
    % Add the next enclosing mask
    nextAdd = find(listToAdd,1);
    imgmask = imgmask | pMaskFcn(pol{nextAdd});
    listToAdd(nextAdd) = false;
    % Subtract any enclosed masks (and then add their children)
    for nextSub = find(tree(:,nextAdd))'
        imgmask = imgmask & ~pMaskFcn(pol{nextSub});
        % Add any children 
        listToAdd(tree(:, nextSub)) = true;
    end
end

end



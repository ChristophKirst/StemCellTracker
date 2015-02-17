function polys = pointsToPolygons(pts, varargin)
%
% shapes = pointsToPolygons(pks, varargin)
%
% description:
%    converts set of points to a set of shape borders using alpha volume
%
% input:
%    pts    array of point coordinates as row vectors
%    param  parameter sturct with entries
%           .radius   radius used in alpha volume detection (100)
%           .dilate   dilate alphavolume with this width ([] = no dilation)
%           .split    split non connected components ([] = false)
%
% output:
%   shapes  cell array of shape borders as arrays of the boundary point coordinates as row vectors
%
% See also: detectAlphaVolume

param = parseParameter(varargin);

% alpha vol
r = getParameter(param, 'radius', 100);
[~, av] = detectAlphaVolume(pts', r);

tri = triangulation(av.tri, pts');

pol = polygonFromTriangulation(tri);

%dilate
di = getParameter(param, 'dilate', []);
if ~isempty(di) && isnumeric(di)
   pol = polygonDiltate(pol, di);
end

%split
sp = getParameter(param, 'split', []);
if ~isempty(sp) && sp
   pol = polygonSplit(pol);
else
   pol = {pol};
end

%convert to polygons
np = length(pol);
polys(np) = ROIPolygon;
for i = 1:np
   polys(i) = ROIPolygon(pol{i});
end

end

% bnd = av.bnd;
% 
% %connected components
% a=sparse(bnd(:,1),bnd(:,2),1);
% a=a+speye(size(a));
% [p,~,r,~]=dmperm(a);
% 
% % sort into groups
% ns = length(r)-1;
% shapes = cell(1,ns);
% for ii=1:ns
%     shapes{ii}=p(r(ii):(r(ii+1)-1));
% end
% shapes(cellfun(@length, shapes)==1)=[];
% ns = length(shapes);
% 
% % split shapes that share a point
% 
% % order the points
% bfrom = bnd(:,1)';
% bto   = bnd(:,2)';
% 
% for s = 1:ns
%    id = shapes{s}(1);
%    
%    pol = zeros(1, length(shapes{s})+1);
%    pol(1) = id;
%    
%    for i = 1:length(shapes{s})
%       eid = find(bfrom == id, 1);
%       id = bto(eid);
%       pol(i+1) = id;
%    end
% 
%    shapes{s} = pts(:,pol);
% end



% % check for holes = shapes inside shapes and remove
% if getParameter(param, 'check', false)
% 
%    shapescheck = true(1, ns);
%    shapeskeep  = true(1, ns);
%    
%    for sc = 1:ns
%       if shapescheck(sc)
%          shc = shapes{sc};
%          for i = find(shapeskeep)
%             if i ~= sc
%                if ispolygonInPolygon(shapes{i}, shc)
%                   shapescheck(i) = 0;
%                   shapeskeep(i)  = 0;
%                end
%             end
%          end
%          %shapescheck(sc) = 0;
%       end
%    end
%    
%    shapes = shapes(shapeskeep);
% end

% dilate
% di = getParameter(param, 'dilate', []);
% if ~isempty(di)
%    for s = 1:length(shapes)
%       shapes{s} = dilatePolygon(shapes{s}, di);
%    end
% end

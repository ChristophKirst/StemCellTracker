function shapes = points2shapes(pts, varargin)
%
% shapes = points2shapes(pks, varargin)
%
% description:
%    converts set of points to a set of shape borders using alpha volume
%
% input:
%    pts    array of point coordinates as row vectors
%    param  parameter sturct with entries
%           .radius   radius used in alpha volume detection (100)
%           .dilate   dilate the shape by this width ([] = no dilation)
%
% output:
%   shapes  cell array of shapeborders as arrays of the boundary point coordinates as row vectors
%
% See also: detectAlphaVolume, dilatePolygon

param = parseParameter(varargin);

% alpha vol
r = getParameter(param, 'radius', 100);
[~, av] = detectAlphaVolume(pts', r);
bnd = av.bnd;

%conneected components
a=sparse(bnd(:,1),bnd(:,2),1);
a=a+speye(size(a));
[p,~,r,~]=dmperm(a);

% sort into groups
ns = length(r)-1;
shapes = cell(1,ns);
for ii=1:ns
    shapes{ii}=p(r(ii):(r(ii+1)-1));
end
shapes(cellfun(@length, shapes)==1)=[];
ns = length(shapes);

% order the points
bfrom = bnd(:,1)';
bto   = bnd(:,2)';

for s = 1:ns
   id = shapes{s}(1);
   
   pol = zeros(1, length(shapes{s})+1);
   pol(1) = id;
   
   for i = 1:length(shapes{s})
      eid = find(bfrom == id, 1);
      id = bto(eid);
      pol(i+1) = id;
   end

   shapes{s} = pts(:,pol);
end

% check for holes = shapes inside shapes and remove

shapescheck = ones(1, ns);
shapeskeep  = ones(1, ns);

for sc = 1:ns
   if shapescheck(sc)
      shc = shapes{sc};
      for i = find(shapeskeep)
         if i ~= sc
            if ispolygonInPolygon(shapes{i}, shc)
               shapescheck(i) = 0;
               shapeskeep(i)  = 0;
            end
         end
      end
      %shapescheck(sc) = 0;
   end
end

shapes = shapes(shapeskeep > 0);

% dilate
bw = getParameter(param, 'dilate', []);
if ~isempty(bw)
   for s = 1:length(shapes)
      shapes{s} = dilatePolygon(shapes{s}, bw);
   end
end

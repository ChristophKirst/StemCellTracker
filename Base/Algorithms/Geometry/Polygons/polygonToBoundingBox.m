function bbox = polygonToBoundingBox(pol)
%
% bbox = polygonToBoundingBox(pol)
%
% description:
%     determines bounding box of polygon 
%
% input: 
%     pol     polygon
%
% output:
%    bbox     rectangular polygon representing the bounding box 
     
if ~iscell(pol) 
   pol = {pol};
end

pmax = max(cell2mat(cellfunc(@(x) max(x, [], 2), pol)), [],2);
pmin = min(cell2mat(cellfunc(@(x) min(x, [], 2), pol)), [],2);

bbox = [pmin, [pmin(1); pmax(2)], pmax, [pmax(1); pmin(2)]];

end
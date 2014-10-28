function pks = stitchPoints(pts, shifts, isizes, varargin)
%
% pks = stitchPoints(peaks, shifts, isizes, varargin)
%
% description:
%       takes points from differnet images together with their shifts and sizes
%       and generates a stitched point set pks
%
% input:
%       pts      list of point coordinates in each image
%       shifts   image shifts
%       isizes   image sizes
%       param    (optional) parameter struct with entries
%                .method    'combine' = join all points together, fast, 
%                           'overwrite'= only add points form a single tile in overlap regions ('combine')

if isempty(pts)
   pks = [];
   return
end

if numel(pts) ~= numel(shifts) 
   error('stitchPeaks: inconsistent peaks and shifts sizes %g ~= %g', numel(pts), numel(shifts));
end

param = getParameter(varargin);

mth = getParameter(param, 'method', 'combine');

% shift peaks
pts = cellfun(@(x,y) x + repmat(y(:), 1, size(x,2)), pts, shifts, 'UniformOutput', false);

switch mth
   case 'combine'
      pks = cell2mat(pts);

   case 'overwrite'
      pks = pts{1};
      for i = 2:numel(pts)
         % clean area
         pks = cleanPoints(pks, shifts{i} + 1, shifts{i} + isizes{i});
         pks = [pks, pts{i}]; %#ok<AGROW>
      end
      
   otherwise
      error('stitchPeaks: no method %s available', mth);
end

end


function pks = cleanPoints(pks, p1, p2) 
   % removes all points in specified area
   p1 = repmat(p1, 1, size(pks,2));
   p2 = repmat(p2, 1, size(pks,2));

   pks(all(and(pks >= p1, pks <= p2), 1)) = [];
end
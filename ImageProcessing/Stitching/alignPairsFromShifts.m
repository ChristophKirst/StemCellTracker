function pairs = alignPairsFromShifts(pairs, shifts, varargin)
%
%  pairs = alignmentPairsFromShifts(pairs, shifts)
%  pairs = alignmentPairsFromShifts(pairs, shifts, nodes)
%
% description:
%    generate AlignmentPair classes consistent with the shifts arranged on a gird 
%

np = length(pairs);

if nargin < 3
   for i = 1:np
      pairs(i).shift = shifts{pairs(i).to} - shifts{pairs(i).from};
   end
else
   % remap nodes to ids in shifts
   nodes = varargin{1};
   [~,remp] = ismember(1:max(nodes), nodes);
   
   for i = 1:np
      pairs(i).shift = shifts{remp(pairs(i).to)} - shifts{remp(pairs(i).from)};
   end
end
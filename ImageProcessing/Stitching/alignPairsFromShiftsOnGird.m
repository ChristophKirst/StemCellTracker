function p = alignPairsFromShiftsOnGird(shifts, varargin)
%
%  p = alignPairsFromShiftsOnGird(shifts)
%
% description:
%    generate AlignmentPair classes consistent with the shifts arragned on a gird 
%

si = size(shifts);
si = padright(si, 3, 1);

np = si(1) * si(2) * (si(3)-1) + si(1) * (si(2)-1) * si(3) + (si(1)-1) * si(2) * si(3);

% if nargin < 2
   p(np) = AlignmentPair();
% else
%    p = varargin{1};
%    if length(p) ~= np
%       error('alignmentPairsFromShiftsOnGird: inconsistent input sizes!');
%    end
% end

i = 1;
for x = 1:si(1)
   for y = 1:si(2)
      for z = 1:si(3)
         if x < si(1)
            p(i).from        = imsub2ind(si, [x,y,z]);
            p(i).to          = imsub2ind(si, [x+1,y,z]);
            p(i).orientation = 1;
            p(i).shift = shifts{p(i).to} - shifts{p(i).from};
            i = i + 1;
         end
         if y < si(2)
            p(i).from = imsub2ind(si, [x,y,z]);
            p(i).to   = imsub2ind(si, [x,y+1,z]);
            p(i).orientation = 2;
            p(i).shift = shifts{p(i).to} - shifts{p(i).from};
            i = i + 1;
         end
         if z < si(3)
            p(i).from = imsub2ind(si, [x,y,z]);
            p(i).to = imsub2ind(si, [x,y,z+1]);
            p(i).orientation = 3;
            p(i).shift = shifts{p(i).to} - shifts{p(i).from};
            i = i + 1;
         end
      end
   end
end


end
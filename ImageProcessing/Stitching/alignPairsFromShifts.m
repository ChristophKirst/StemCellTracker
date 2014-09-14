function pairs = alignPairsFromShifts(pairs, shifts)
%
%  pairs = alignmentPairsFromShifts(pairs, shifts)
%
% description:
%    generate AlignmentPair classes consistent with the shifts arragned on a gird 
%

np = length(pairs);
for i = 1:np
   pairs(i).shift = shifts{pairs(i).to} - shifts{pairs(i).from};
end

end
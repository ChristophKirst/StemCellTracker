function pairs = alignPairsFromOverlap(shifts, isizes)
%
% pairs = overlappingPairs(shifts, isizes)
%
% description:
%     calculate alignment pairs from overlapping images specified by shifts and image sizes
%
% input:
%    shifts   shifts of theimages
%    isizes   sizes of the images
%
% output:
%    pairs    array of Alignment pairs 


n = numel(shifts);

pairs = [];
for i = 1:n
   sh = shifts{i};
   si = isizes{i};
   
   for j = i+1:n
      % check overlap 
      sh2 = shifts{j};
      si2 = isizes{j};
      if all(max(sh,sh2) < min(sh+si, sh2+si2))
         p = AlignmentPair;
         p.from = i;
         p.to   = j;
         p.shift = sh2 - sh;
         pairs = [pairs, p]; %#ok<AGROW>
      end
   end
end
   
   
end

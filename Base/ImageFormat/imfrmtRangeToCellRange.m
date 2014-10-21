function range = imfrmtRangeToCellRange(range)
%
% tgir = imfrmtRangeToCellRange(refTagRange, tagRange)
% 
% description:
%    converts all numerical arrays of indices in range to cells

tnames = fieldnames(range);

for i = 1:length(tnames)
   f = tnames{i}; 

   v = range.(f);
   if iscell(v)
      range.(f) = v;
   else
      range.(f) = num2cell(v);
   end
end

end
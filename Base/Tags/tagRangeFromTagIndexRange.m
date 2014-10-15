function tagRange = tagRangeFromTagIndexRange(refTagRange, iRange)
%
% tgir = tagIndexRangeToTagRange(refTagRange, tagRange)
% 
% description:
%    converts a index range iRange to a tag range
%    numeric arrays are converted to cells
%    numeric arrays that have a cellstr reference in refTagRange are converted to a cellstr using the names


tagRange = iRange;

tnames = fieldnames(iRange);

for i = 1:length(tnames)
   f = tnames{i};
   
   if isfield(refTagRange, f) 
      vref = refTagRange.(f);
      if iscellstr(vref)
         v = tagRange.(f);
         tagRange.(f) = vref(v);
      else
         tagRange.(f) = num2cell(tagRange.(f));
      end
   else
      tagRange.(f) = num2cell(tagRange.(f));
   end
end

end
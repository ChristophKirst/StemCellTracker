function tagIndexRange = tagRangeToTagIndexRange(refTagRange, tagRange)
%
% tgir = tagRangeToTagIndexRange(refTagRange, tagRange)
% 
% description:
%    converts a tag range tagRange to a range compatible with the imfrmt range
%    cells of numeric indices are converted to arrays
%    cells of chars are converted to arrays of indices using the reference tag range refTagRange


tagIndexRange = tagRange;

tnames = fieldnames(tagRange);

for i = 1:length(tnames)
   f = tnames{i};
   
   v = tagRange.(f);
   if iscellstr(v)
      if isfield(refTagRange, f)
         vref = refTagRange.(f);
         if ~iscellstr(vref)
            error('tagRangeToTagIndexRange: type mismatch in field %s in reference and tan range!', f);
         else
            tagIndexRange.(f) = find(ismember(vref, v));
         end
      else
         error('tagRangeToTagIndexRange: cannot find field %s in reference range!', f);
      end
   else
      if iscell(v)
         tagIndexRange.(f) = cell2mat(v);
      else
         tagIndexRange.(f) = v;
      end
   end
end

end
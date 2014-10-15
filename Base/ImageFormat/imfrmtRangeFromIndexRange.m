function range = imfrmtRangeFromIndexRange(refRange, range)
%
% tgir = imfrmtRangeFromIndexRange(refTagRange, tagRange)
% 
% description:
%    converts a index range to a tag range
%    numeric arrays are converted to cells
%    numeric arrays that have a cellstr reference in refTagRange are converted to a cellstr using the names

tnames = fieldnames(range);

refnames = fieldnames(refRange);
refnamesl = lower(refnames);

for i = 1:length(tnames)
   f = tnames{i}; 
   fl = lower(f);
   
   id = find(ismember(refnamesl, fl), 1, 'last');
   if ~isempty(id)
      rf = refnames{id};
      vref = refRange.(rf);
      if ischar(vref)
         vref = {vref};
      end
      v = range.(f);

      if f ~= rf
         range.(f) = vref(end - v + 1);
      else
         range.(f) = vref(v);
      end

   else
      v = range.(f); 
      if iscell(v)
         range.(f) = v;
      else
         range.(f) = num2cell(range.(f));
      end
   end
end

end
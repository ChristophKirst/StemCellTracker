function range = imfrmtIndexRangeFromIndexRange(refRange, range)
%
% tgir = imfrmtIndexRangeFromIndexRange(refTagRange, tagRange)
% 
% description:
%    converts a index range to a tag range
%    numeric arrays / cells are converted to arrays
%    no cellstrs are allowed

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
%       if ischar(vref)
%          vref = {vref};
%       end
      v = range.(f);

      if f ~= rf
         range = rmfield(range, f);
         range.(rf) = vref(end - v + 1);
      else
         range.(f) = vref(v);
      end

   else
      v = range.(f); 
      if iscell(v)
         range.(f) = cell2mat(v);
      end
   end
end

end
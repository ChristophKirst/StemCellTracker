function range = imfrmtRangeFromIndexRange(refRange, range)
%
% tgir = imfrmtRangeFromIndexRange(refTagRange, tagRange)
% 
% description:
%    converts a index range to a tag range
%    numeric arrays that have a cellstr reference in refTagRange are converted to a cellstr

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
         range = rmfield(range, f);
         range.(rf) = vref(end - v + 1);
      else
         range.(f) = vref(v);
      end

   else
%       v = range.(f); 
%       if iscell(v)
%          range.(f) = v;
%       else
%          range.(f) = num2cell(v);
%       end
   end
end

end
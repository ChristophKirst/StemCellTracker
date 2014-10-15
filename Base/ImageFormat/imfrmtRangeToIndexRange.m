function range = imfrmtRangeToIndexRange(refRange, range)
%
% tgir = imfrmtRangeToIndexRange(refRange, range)
% 
% description:
%    converts a coordinate range to a range compatible with the imfrmt range
%    cells of numeric indices are converted to arrays
%    entries present in refRange are converted to arrays indices
%    changes due to flips in the format are taken care of

tnames = fieldnames(range);

refnames = fieldnames(refRange);
refnamesl = lower(refnames);

for i = 1:length(tnames)
   f = tnames{i};
   fl = lower(f);
 
   v = range.(f);
   if ischar(v)
      v = {v};
   end
   
   id = find(ismember(refnamesl, fl), 1, 'last');
   if ~isempty(id)
      rf = refnames{id};

      vref = refRange.(rf);
      if ischar(vref)
         vref = {vref};
      end

      if iscellstr(v)
         if iscellstr(vref)
            v = find(ismember(vref, v));            
         else
            error('imfrmtRangeToIndexRange: type mismatch in field %s in reference and range!', f);
         end
      else
         if iscell(v)
            v = cell2mat(v);
         end         
         if isnumeric(vref) && isnumeric(v)
            if iscell(vref)
               vref = cell2mat(vref);
            end
            
            [id, v] = ismember(v, vref); % indices should all be covered
            if any(~id)
               error('imfrmtRangeToIndexRange: values in field %s out of reference range!', f);
            end
         elseif ~(iscellstr(vref) && isnumeric(v))
            % if vref is cellstr but v is numeric then just take v
            error('imfrmtRangeToIndexRange: type mismatch in field %s in reference and range!', f);
         end
      end
  

      if rf ~= f
         range.(f) = length(vref) - v + 1;
      else
         range.(f) = v;
      end

   else % not in ref range
      if iscell(v)
         range.(f) = cell2mat(v);
      else
         range.(f) = v;
      end
   end
end

end
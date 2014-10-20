function range = imfrmtRangeToIndexRange(refRange, varargin)
%
% range = imfrmtRangeToIndexRange(refRange, range)
% 
% description:
%    converts a coordinate range to a range compatible with the imfrmt range
%    cells of numeric indices are converted to arrays
%    entries present in refRange are converted to array indices
%    changes due to flips in the format produce error 

range = parseParameter(varargin);

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

%       v
%       vref
      if iscellstr(v)
         if iscellstr(vref)
            vi = ismember(vref, v);            
         else
            error('imfrmtRangeToIndexRange: type mismatch in field %s in reference and range!', f);
         end
         if sum(vi) ~= length(v)
            error('imfrmtRangeToIndexRange: %s keys %s do not exists in range %s!', rf, var2char(v(~vi)), var2char(vref));
         else
            v = find(vi);
         end
         
      else
         if iscell(v)
            v = cell2mat(v);
         end
         if iscell(vref)
            vref = cell2mat(vref);
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
            % if vref is cellstr but v is numeric -> check if in range
            if any(v < 1) || any(v > length(vref))
               error('imfrmtRangeToIndexRange: index out of bounds in field %s!', f);
            end
         end
      end
  

      if rf ~= f
         %range = rmfield(range, f);
         %length(vref)
         %range.(rf) = length(vref) - v + 1;
         error('imfrmtRangeToIndexRange: inversion axis not recommended for non index ranges!');
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
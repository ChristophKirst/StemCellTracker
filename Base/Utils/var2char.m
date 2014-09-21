function c = var2char(cs)
%
% c = var2char(cs)
%
% description: converts a cell array of strings to a char for printing to screen

% note: as matlab, here columns are used as primary index

if ischar(cs)
   c = ['''' , cs, ''''];
   return
elseif ~iscell(cs) && isscalar(cs)
   if isnumeric(cs) || islogical(cs)
      c = num2str(cs);
   else
      c = ['<' class(cs) , '>'];
   end
else
   if isempty(cs)
      if iscell(cs)
         c = '{}';
      else
         c ='[]';
      end
      return
   end
 
   si = size(cs);
   di = length(si);
   if di == 2 && si(2) == 1
      di = 1;
      si(2) = [];
   end
   
   cc = {};
   cdi = 0;
   for d = 1:di-1
      cc = [cc, {1:si(d)}]; %#ok<AGROW>
      if si(d) > 1
         cdi = 1;
      end
   end
   
   if iscell(cs)
      if isempty(cs)
         c ='{}';
         return
      end
      
      c = '{';
      if cdi
         for i = 1:si(end)-1
            c = [c var2char(cs(cc{:}, i)) ', ']; %#ok<AGROW>
         end
         c = [c var2char(cs(cc{:}, si(end))) '}'];
      else
         for i = 1:si(end)-1
            c = [c var2char(cs{cc{:}, i}) ', ']; %#ok<AGROW>
         end
         c = [c var2char(cs{cc{:}, si(end)}) '}'];
      end 
   else
      c = '[';
      for i = 1:si(end)-1
         c = [c var2char(cs(cc{:}, i)) ', ']; %#ok<AGROW>
      end
      c = [c var2char(cs(cc{:}, si(end))) ']'];
   end
end
   
   


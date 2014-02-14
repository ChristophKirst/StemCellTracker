function b = isimarishead(id)
%
% b = isimarishead(id)
%
% description:
%    checks if id is Imaris application, an integer i or a string that is 
%    convertible to a Imaris application id
%
% input:
%    id    id to check
%
% output:
%    b     true / false
%
% See also: isimarisid, isimaris

b = 0;

if isa(id, 'Imaris.IApplicationPrxHelper')
   b = 1;
   return
end

if isnumeric(id) 
   if isscalar(id)
      b = 1;
   end
   return
end

if ischar(id) && ~isnan(str2double(id))
   b = 1;
end

end
      
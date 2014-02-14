function b = isimaris(object)
%
% b = isimaris(object)
%
% description:
%   checks if b is an Imaris application object
%
% input:
%   object      object to check
%
% output:
%   b           true / false
%
% See also: isimarisid, imarisinitialize, imarisstart

b = isa(object,  'Imaris.IApplicationPrxHelper');

end
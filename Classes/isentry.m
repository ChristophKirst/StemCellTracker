function ie = isentry(obj, name)
%
% bool = isentry(obj, name)
%
% description:
%   checks if name is a field, property or method of obj
%
% input:
%   obj     object
%   name    field, property or method name of obj
%

if length(obj) > 1
   ie = isfield(obj, name) || isprop(obj(1), name) || ismethod(obj, name);
else
   ie = isfield(obj, name) || isprop(obj, name) || ismethod(obj, name);
end
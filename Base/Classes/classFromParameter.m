function obj = classFromParameter(obj, prefix, varargin)
%
% obj = classFromParameter(obj, prefix, varargin)
%
% description:
%     initializes class from parameter using a prefrix

param = parseParameter(varargin);

pnames = fieldnames(param);

for i = 1:length(pnames)
   if isprop(obj, [prefix, pnames{i}])
      obj.([prefix, pnames{i}]) = param.(pnames{i});
   else
      warning('%s: unknown property name: %s ', class(obj), pnames{i})
   end
end

end
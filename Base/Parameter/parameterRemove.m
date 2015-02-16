function param = parameterRemove(param, fnames)
%
% param = parameterRemove(param, fnames)
%
%
% description:
%    removes all fields fnames
%
% input:
%    param   parameter struct
%    fnames  field names to remove
%
% output:
%    param    parmeter struct with fnames removed

pnames = fieldnames(param);
rmnames = intersect(pnames, fnames);
param = rmfield(param, rmnames);

end
function param = parameterFilter(param, fnames)
%
% param = parameterFilter(param, fnames)
%
%
% description:
%    removes all fields not in fnames
%
% input:
%    param   parameter struct
%    fnames  field names to keep
%
% output:
%    param    parmeter struct with at most fnames fields

pnames = fieldnames(param);
rmnames = setdiff(pnames, fnames);
param = rmfield(param, rmnames);

end
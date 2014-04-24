function b = isemptystruct(s)

b = isempty(s) || (isstruct(s) && isempty(fieldnames(s)));

end
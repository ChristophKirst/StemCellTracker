function b = isemptystruct(s)

b = isempty(s) || isempty(fieldnames(s));

end
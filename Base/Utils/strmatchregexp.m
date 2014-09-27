function s = strmatchregexp(expr, str)
%
% s = strmatchregexp(expr, str)
%
% description:
%    checks if regular experssion matches one of the srings in str
%

if ~iscell(str)
   str = {str};
end

s = ~cellfun(@isempty, regexp(expr, str));

end
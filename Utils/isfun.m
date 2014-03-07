function b = isfun(f)
%
% b = isfun(f)
%
% description:
%    returns ture if f is a function handle
%
% input:
%    f  variable to test
%
% output:
%    b  true if f is function_handle

b = isa(f, 'function_handle');

end
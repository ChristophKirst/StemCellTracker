function b = isfile(fname)
%
% b = isfile(fname)
%

%note: exists may be to smart
b = (exist(fname, 'file') == 2);

end
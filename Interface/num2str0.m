function str = num2str0(num, n, pad)
%
% str = num2str0(num, n, pad)
%
% decription: 
%    converts num to a string of length n padding form the left with pad
%
% input:
%    num   integer to convert
%    n     number of digits
%    pad   (optional) padding character ('0')

if nargin < 3
   pad = '0';
end

str = num2str(num);
l =length(str);

if l > n
   error('num2str0: number string %s is langer than %g', str, n)
end

str =[repmat(pad, 1, n - l), str];

end




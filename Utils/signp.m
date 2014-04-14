function sgn = signp(val)
%
% sgn = signp(val)
%
% description:
%    returns sign(val) but with result 0 replaced by 1

sgn = sign(val);
sgn(sgn==0) = 1;

end
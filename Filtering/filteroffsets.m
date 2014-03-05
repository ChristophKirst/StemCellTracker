function [ol, or] = filteroffsets(ksize)
%
% [ol, or] = filteroffsets(hsize)
%
% description:
%    returns bw image with pixels at points set to white
%
% input:
%    ksize     input size of filter
%
% output:
%    ol,or     radial offsets
%
% note:
%    for even numbers the filer is moved towards the left in accordance wirh matlab routines
%
% See also: fspecial, fspecial3


for i=length(ksize):-1:1   
   if mod(ksize(i),2) == 0
      ol(i) = ksize(i)/2;
      or(i) = ksize(i)/2-1;
   else
      ol(i) = (ksize(i)-1)/2;
      or(i) = ol(i);
   end
end

end
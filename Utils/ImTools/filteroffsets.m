function [ol, or] = filteroffsets(hsize)
%
% [ol, or] = filteroffsets(hsize)
%
% description:
%    returns bw image with pixels at points set to white
%
% input:
%    hsize     input size of filter
%
% output:
%    ol,or     radiaul offsets
%
% note:
%    for even numbers matlab moves the filer towards the left
%
% See also: fspecial, fspecial3


for i=length(hsize):-1:1   
   if mod(hsize(i),2) == 0
      ol(i) = hsize(i)/2;
      or(i) = hsize(i)/2-1;
   else
      ol(i) = (hsize(i)-1)/2;
      or(i) = ol(i);
   end
end

end
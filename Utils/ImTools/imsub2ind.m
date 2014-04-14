function idx = imsub2ind(siz, sub)
%
% out = imsub2ind(siz, sub)
%
% description: 
%     converst array of position coordinates to the linear indices
%     assuming image dimensions d
%
% input:
%     siz      size vector
%     sub      indices as array of the form [idim1, idim2, ...] (e.g. as returned by imind2sub)
%
% output:
%     idx      linear indices
%
% See also: sub2ind, imind2sub


siz = double(siz);
lensiz = length(siz);
s = size(sub,2); 

if s ~= lensiz
   error('imsub2ind: dimension mimatch');
end

%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
idx = 1;

for i = 1:s
   v = sub(:,i);
   if (any(v(:) < 1)) || (any(v(:) > siz(i)))
      %Verify subscripts are within range
      error('imsub2ind: index out of range');
   end
   idx = idx + (v-1)*k(i);
end
    
end




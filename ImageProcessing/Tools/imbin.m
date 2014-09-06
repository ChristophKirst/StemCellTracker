function imgb = imbin(img, binsize, meth)
%
% imgb = imbin(img, binsize, meth)
%
% description:
%    bins the pixels using binsize and then calcualtes the new pixel value using the meth: mean or median


if nargin < 3
   meth = @mean;
end

dim = ndims(img);

if length(binsize) == 1 
   binsize = ones(1,dim) * binsize;
end
%binsize
   
isize = size(img);
isizebin = ceil(isize ./ binsize - eps);
imgb = zeros(isizebin);

if dim == 2
   for i = 1:isizebin(1)
      for j = 1:isizebin(2)
         dat = img((i-1)*binsize(1)+1:min(i*binsize(1), isize(1)), (j-1)*binsize(2)+1:min(j*binsize(2), isize(2)));
         imgb(i,j) = meth(dat(:));
      end
   end
else
   for i = 1:isizebin(1)
      for j = 1:isizebin(2)
         for k = 1:isizebin(3)
            dat = img((i-1)*binsize(1)+1:min(i*binsize(1), isize(1)), (j-1)*binsize(2)+1:min(j*binsize(2), isize(2)),...
                      (k-1)*binsize(3)+1:min(k*binsize(3), isize(3)));
            imgb(i,j,k) = meth(dat(:));
         end
      end
   end
end

end
         




% [m,n]=size(M); %M is the original matrix
% 
% M=sum( reshape(M,p,[]) ,1 );
% M=reshape(M,m/p,[]).'; %Note transpose
% 
% M=sum( reshape(M,q,[]) ,1),
% M=reshape(M,n/q,[]).'; %Note transpose 



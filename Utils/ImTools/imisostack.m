function isostack = imisostack(stack, zslices, method)
%
% stack = imisostack(image,  zslices)
%
% description:
%     returns isotropic stack by copying z-slices to match xyspacing
%
% input:
%     stack    image stack
%     zslices  number of zslices to insert between two subsequent slices
%     method   (optional) mehod to detemine intermediate images ('copy' (copy image), 
%                         'liner' (linaer interpol), 'mean', or fucntion) 
%                         'max' (maximal intensity), 'min' (minimal intensity) ('copy')
%
% output:
%     isostack  new isotropic stack

siz = size(stack);

if length(siz) ~=3
   error('imisostack: expect 3d stack of gray scale images.')
end

if nargin < 3
   method = 'copy';
end

if zslices > 0
   zslices = zslices + 1; %zslices are now total number of images between z and z+1 in original image
   if ischar(method)
      switch method
         case 'copy'
            isostack = zeros(siz(1), siz(2), zslices * siz(3));
            for i = 1:siz(3);
               isostack(:,:,(i-1)*zslices+1:i*zslices) = repmat(stack(:,:,i),[ 1,1,zslices]);
            end        
  
         case 'linear'
            isostack = zeros(siz(1), siz(2), zslices * (siz(3)-1) + 1);
            for i = 1:siz(3)-1;
               for z = 1:zslices
                  zz = (i-1)*zsclices + z;
                  isostack(:,:,zz) = (zslices-z+1)/zslices * double(stack(:,:,i)) + (z-1)/zslices * double(stack(:,:,i+1));
               end 
            end
            isostack(:,:,end) = stack(:,:,end);
  
         case 'mean'
            isostack = zeros(siz(1), siz(2), zslices * (siz(3)-1) + 1);
            for i = 1:siz(3)-1;
               immean = (double(stack(:,:,i)) + double(stack(:,:,i+1)))/2;
               isostack(:,:,(i-1)*zslices+1) = stack(:,:,i);
               isostack(:,:,(i-1)*zslices+2:i*zslices) = repmat(immean,[1,1,zslices-1]); 
            end
            isostack(:,:,end) = stack(:,:,end);
            
         case 'max'
            isostack = zeros(siz(1), siz(2), zslices * (siz(3)-1) + 1);
            for i = 1:siz(3)-1;
               immax = min(stack(:,:,i), stack(:,:,i+1));
               isostack(:,:,(i-1)*zslices+1) = stack(:,:,i);
               isostack(:,:,(i-1)*zslices+2:i*zslices) = repmat(immax,[1,1,zslices-1]); 
            end
            isostack(:,:,end) = stack(:,:,end);
            
         case 'min'
            isostack = zeros(siz(1), siz(2), zslices * (siz(3)-1) + 1);
            for i = 1:siz(3)-1;
               immin = min(stack(:,:,i), stack(:,:,i+1));
               isostack(:,:,(i-1)*zslices+1) = stack(:,:,i);
               isostack(:,:,(i-1)*zslices+2:i*zslices) = repmat(immin,[1,1,zslices-1]); 
            end
            isostack(:,:,end) = stack(:,:,end);
                 
       otherwise
         error('imisostack: not a valid method: %s', method);
   
      end
   elseif isfun(method)            
      isostack = zeros(siz(1), siz(2), zslices * (siz(3)-1) + 1);
      for i = 1:siz(3)-1;
         for z = 1:zslices
            zz = (i-1)*zsclices + z;
            isostack(:,:,zz) = method(double(stack(:,:,i)), double(stack(:,:,i+1)), z, zslices);
         end
      end
      isostack(:,:,end) = stack(:,:,end);
   else
      error('imisostack: method not valid!')
   end

elseif zslices < 0 % remove zslices
   zslices = floor(abs(zslices) + 1);

   if ischar(method) && strcmp(method, 'mean')
      if rem(siz(3)+1, zslices) ~= 0
         error('imisostack: number of z slices should be multiple of reduction factor!')
      end
      isostack = zeros(siz(1), siz(2),  (siz(3)+1)/ zslices);
      
   else
      if rem(siz(3), zslices) ~= 0
         error('imisostack: number of z slices should be multiple of reduction factor!')
      end
      
      isostack = zeros(siz(1), siz(2),  siz(3)/ zslices);
   end

   if ischar(method) && strcmp(method, 'mean')
      for i = 1:((siz(3)+1)/zslices - 1)
         isostack(:,:,i) = mean(stack(:,:,(i-1)*zslices + 1:i*zslices),3);
      end
      isostack(:,:,end) = stack(:,:,end);
   else   % subsampling
      for i = 1:((siz(3))/zslices)
         isostack(:,:,i) = stack(:,:,(i-1)*zslices+1);
      end
   end
      
end


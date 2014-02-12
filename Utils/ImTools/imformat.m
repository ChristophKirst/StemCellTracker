function format = imformat(data)
%
% format = imformat(data)
%
% description:
%     determines the format of the data
%     color dimensions are detected as having size 3
% 
% input:
%     data    image data
%
% output:
%     format:
%     'hw'   = h x w matrix = 'hw' (2D grayscale)
%     'hwc'  = h x w x c matrix (2D rgb )
%     'hwl'  = h x w x l matrix (3D)
%     'hwcl' = h x w x c x l matrix (3D color image, matlab ordering)
%     'hwlc' = h x w x l x c matrix (3D color image)
%     'hwlt'  = h x w x l x t matrix (4D)
%     'hwclt' = h x w x c x l x t matrix (4D color image, matlab ordering)
%     'hwlct' = h x w x l x c x t matrix (4D color image)
%     ''  = not supported format
%     h = height, w = width, l = length, c = color, t = time
%
% note:
%     'hwl' = 'hwt', 'hwcl' = 'hwct', 'hwlc' = 'hwtc'


dim = ndims(data);

format = '';

if dim < 2 || dim > 4
   return
end

switch dim
   case 2
      format = 'hw';
   case 3
      if size(data,3) == 3
         format = 'hwc';
      else
         format = 'hwl';
      end
   case 4
      if size(data,3) == 3
         format = 'hwcl';
      elseif size(data,4) == 3
         format = 'hwlc';
      else
         format = 'hwlt';
      end
   case 5
      if size(data,3) == 3
         format = 'hwclt';
      elseif size(data,4) == 3
         format = 'hwlct';
      elseif size(data,5) == 3
         format = 'hwltc';
      end
end

end
      
   





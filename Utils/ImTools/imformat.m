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
%     'hwc'  = h x w x c matrix (2D multi channel )
%     'hwl'  = h x w x l matrix (3D)
%     'hwcl' = h x w x c x l matrix (3D multi channel image, matlab ordering)
%     'hwlc' = h x w x l x c matrix (3D multi channel image)
%     'hwlt'  = h x w x l x t matrix (4D grayscale)
%     'hwclt' = h x w x c x l x t matrix (4D multi channel image, matlab ordering)
%     'hwlct' = h x w x l x c x t matrix (4D multi channel image, time last)
%     'hwlct' = h x w x l x t x c matrix (4D multi channel image, channel last)
%     ''  = not supported format
%     h = height, w = width, l = length, c = color, t = time
%
% note:
%     'hwl' = 'hwt', 'hwcl' = 'hwct', 'hwlc' = 'hwtc'

cmaxsize = 4;  % max size for channel dimension

dim = ndims(data);

format = '';

if dim < 2 || dim > 5
   return
end

switch dim
   case 2
      format = 'hw';
   case 3
      if size(data,3) <= cmaxsize
         format = 'hwc';
      else
         format = 'hwl';
      end
   case 4
      if size(data,3) <= cmaxsize
         format = 'hwcl';
      elseif size(data,4) <= cmaxsize
         format = 'hwlc';
      else
         format = 'hwlt';
      end
   case 5
      if size(data,3) <= cmaxsize
         format = 'hwclt';
      elseif size(data,4) <= cmaxsize
         format = 'hwlct';
      elseif size(data,5) <= cmaxsize
         format = 'hwltc';
      else
         format = 'hwlct';
      end     
end

end
      
   





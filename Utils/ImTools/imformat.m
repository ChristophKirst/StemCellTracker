function format = imformat(data)
%
% format = imformat(data)
%
% description:
%     determines the format of the data
%     channel dimension is detected as dimensions having size <= 5
% 
% input:
%     data    image data
%
% output:
%     format:
%     'pq'   = p x q matrix = 'pq' (2D grayscale)
%     'pqwc'  = p x q x c matrix (2D multi channel )
%     'pql'  = p x q x l matrix (3D)
%     'pqcl' = p x q x c x l matrix (3D multi channel image, matlab ordering)
%     'pqlc' = p x q x l x c matrix (3D multi channel image)
%     'pqlt'  = p x q x l x t matrix (4D grayscale)
%     'pqclt' = p x q x c x l x t matrix (4D multi channel image, matlab ordering)
%     'pqlct' = p x q x l x c x t matrix (4D multi channel image, time last)
%     'pqlct' = p x q x l x t x c matrix (4D multi channel image, channel last)
%     ''  = not supported format
%     p = x pixel coordinate, q = z pixel coordinate, l = z pixel coordinate, c = color, t = time
%
% note:
%     'pql' = 'pqt', 'pqcl' = 'pqct', 'pqlc' = 'pqtc'

cmaxsize = 4;  % max size for channel dimension

dim = ndims(data);

format = '';

if dim < 2 || dim > 5
   return
end

switch dim
   case 2
      format = 'pq';
   case 3
      if size(data,3) <= cmaxsize
         format = 'pqc';
      else
         format = 'pql';
      end
   case 4
      if size(data,3) <= cmaxsize
         format = 'pqcl';
      elseif size(data,4) <= cmaxsize
         format = 'pqlc';
      else
         format = 'pqlt';
      end
   case 5
      if size(data,3) <= cmaxsize
         format = 'pqclt';
      elseif size(data,4) <= cmaxsize
         format = 'pqlct';
      elseif size(data,5) <= cmaxsize
         format = 'pqltc';
      else
         format = 'pqlct';
      end     
end

end
      
   





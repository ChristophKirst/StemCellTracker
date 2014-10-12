function frmt = imfrmtFormatFromSize(isize, cmaxsize)
%
% format = imfrmtFormatFromSize(isize)
%
% description:
%     tries to guess the format of an image with size isize
%     channel dimension is detected as having dimensions <= cmaxsize = 4
% 
% input:
%     isize    image size
%     cmaxsize (optional) maximal channel for colors (4)
% 
% output:
%     frmt     guessed format, '' if failed to determine the format
%
% note:
%     'XYZ' = 'XYT', 'XYCZ' = 'XYCT', etc, Z is given preference over T
%
% See also: imfrmtFormat

if nargin < 2
   cmaxsize = 3;  % max size for channel dimension
end

dim = length(isize);

frmt = ''; 

if dim < 2 || dim > 5
   return
end

switch dim
   case 2
      frmt = 'XY';
   case 3
      if isize(3) <= cmaxsize
         frmt = 'XYC';
      else
         frmt = 'XYZ';
      end
   case 4
      if isize(3) <= cmaxsize
         frmt = 'XYCZ';
      elseif isize(4) <= cmaxsize
         frmt = 'XYZC';
      else
         frmt = 'XYZT';
      end
   case 5
      if isize(3) <= cmaxsize
         frmt = 'XYCZT';
      elseif isize(4) <= cmaxsize
         frmt = 'XYZCT';
      elseif isize(5) <= cmaxsize
         frmt = 'XYZTC';
      else
         frmt = 'XYZCT';
      end     
end

end
      
   





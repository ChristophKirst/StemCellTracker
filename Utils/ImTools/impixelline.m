function img = impixelline(img, pos1, pos2, val)
%
% imgln = impixelline(img, pos1, pos2, val)
%
% description:
%      changing the pixel values to val along a line between pos1 and pos2
%
% input:
%      img          non-color 2d image
%      pos1, pos2   points to connect.
%      val          value of pixel along line
% 
% output:
%      imgln        image with added line
%
% todo: 
%      extend to 3d images

x0 = pos1(1); y0 = pos1(2);
x1 = pos2(1); y1 = pos2(2);

isize = size(img);
%if x0 > isize(1) || x1 > isize(1) || y0 > isize(2) || y1 > isize(2)
%   error('impixelline: positions larger than image size');
%end
   

if checkxy(x0,y0,isize); img(x0, y0) = val; end
if checkxy(x1,y1,isize); img(x1, y1) = val; end

if abs(x1 - x0) <= abs(y1 - y0)
   if y1 < y0
      k = x1; x1 = x0; x0 = k;
      k = y1; y1 = y0; y0 = k;
   end
   if (x1 >= x0) && (y1 >= y0)    
      dy = y1-y0; dx = x1-x0;
      p = 2*dx; n = 2*dy - 2*dx; tn = dy;
      while (y0 < y1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; x0 = x0 + 1;
         end
         y0 = y0 + 1;        
         %[x0, y0] = checkxy(x0,y0, isize); 
         if checkxy(x0,y0,isize); img(x0, y0) = val; end
      end
   else
      dy = y1 - y0; dx = x1 - x0;
      p = -2*dx; n = 2*dy + 2*dx; tn = dy;
      while (y0 <= y1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; x0 = x0 - 1;
         end
         y0 = y0 + 1;
         %[x0, y0] = checkxy(x0,y0, isize); 
         if checkxy(x0,y0,isize); img(x0, y0) = val; end
      end
   end
else
   if x1 < x0
      k = x1; x1 = x0; x0 = k;
      k = y1; y1 = y0; y0 = k;
   end
   if (x1 >= x0) && (y1 >= y0)
      dy = y1 - y0; dx = x1 - x0;
      p = 2*dy; n = 2*dx-2*dy; tn = dx;
      while (x0 < x1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; y0 = y0 + 1;
         end
         x0 = x0 + 1; 
         %[x0, y0] = checkxy(x0,y0, isize); 
         if checkxy(x0,y0,isize); img(x0, y0) = val; end
      end
   else
      dy = y1 - y0; dx = x1 - x0;
      p = -2*dy; n = 2*dy + 2*dx; tn = dx;
      while (x0 < x1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; y0 = y0 - 1;
         end
         x0 = x0 + 1;
         %[x0, y0] = checkxy(x0,y0, isize); 
         if checkxy(x0,y0,isize); img(x0, y0) = val; end
      end
   end
end

img = imextract(img, [1, 1], isize);

end



function chk = checkxy(x0, y0, isize)
   chk = false;
   if (x0 > isize(1) || x0 < 1)
      return
   end
   
   if (y0 > isize(2) || y0 < 1)
      return
   end
   
   chk = true;
end
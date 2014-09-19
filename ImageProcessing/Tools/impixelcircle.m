function img = impixelcircle(img, cent, radius, theta1, theta2, val)
%
% imgln = impixelcircle(img, rcent, theta1, theta2, val)
%
% description:
%      changing the pixel values to val along a line between pos1 and pos2
%
% input:
%      img            non-color 2d image
%      cent, radius   center and radius of circle
%      theta1, theta2 circle segment 
%      val            value of pixel along line
% 
% output:
%      img            image with added circle
%

isize = size(img);

xc = int16(cent(1));
yc = int16(cent(2));

x = int16(0);
y = int16(radius);
d = int16(1 - radius);

if check(xc, yc+y, xc, yc, theta1, theta2, isize); img(xc, yc+y) = val; end
if check(xc, yc-y, xc, yc, theta1, theta2, isize); img(xc, yc-y) = val; end
if check(xc+y, yc, xc, yc, theta1, theta2, isize); img(xc+y, yc) = val; end
if check(xc-y, yc, xc, yc, theta1, theta2, isize); img(xc-y, yc) = val; end

while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    if check( x+xc,  y+yc, xc, yc, theta1, theta2, isize); img( x+xc,  y+yc) = val; end
    if check( y+xc,  x+yc, xc, yc, theta1, theta2, isize); img( y+xc,  x+yc) = val; end
    if check( y+xc, -x+yc, xc, yc, theta1, theta2, isize); img( y+xc, -x+yc) = val; end
    if check( x+xc, -y+yc, xc, yc, theta1, theta2, isize); img( x+xc, -y+yc) = val; end
    if check(-x+xc, -y+yc, xc, yc, theta1, theta2, isize); img(-x+xc, -y+yc) = val; end
    if check(-y+xc, -x+yc, xc, yc, theta1, theta2, isize); img(-y+xc, -x+yc) = val; end
    if check(-y+xc,  x+yc, xc, yc, theta1, theta2, isize); img(-y+xc,  x+yc) = val; end
    if check(-x+xc,  y+yc, xc, yc, theta1, theta2, isize); img(-x+xc,  y+yc) = val; end
end

end



function chk = check(x0, y0, xc, yc, th1, th2, isize)
   chk = false;
   if (x0 > isize(1) || x0 < 1)
      return
   end

   if (y0 > isize(2) || y0 < 1) 
      return
   end
   
   th = atan2(double(y0-yc),double(x0-xc));  
   if th < 0
      th = th + 2 * pi;
   end
   
   thr = th - th1;
   if thr < 0
      thr = thr + 2 * pi;
   end
   
   th21 = th2 - th1;
   
   %{x0, y0, th, thr, th21, thr > th21}
   
   if (thr > th21)
      return
   end
   
   chk = true;
end
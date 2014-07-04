function a = overlapDiskSegmentTriangle(diskCenter, diskRadius, diskTheta1, diskTheta2, triangleCorner1, triangleCorner2, triangleCorner3)
% 
% a = overlapDiskTriangle(diskCenter, diskRadius,  diskTheta1, diskTheta2, triangleCorner1, triangleCorner2, triangleCorner3)
% 
% description:
%     calculates the overlap between a disk and a triangle
%
% input:
%    diskCenter, diskRadius     center and radius of disk
%    diskTheta*                 start and end angle of the disk segment in [0, 2 pi], 0 is x axis, theta is measured counter clockwise
%    triangleCorner*            corners of the triangle
%
% output:
%    a                         overlap


% relative coords

x1 = triangleCorner1 - diskCenter;
x2 = triangleCorner2 - diskCenter;
x3 = triangleCorner3 - diskCenter;


% overlaps for all 

[x, yl, yr, t1, t2] = standardCoords(x1, x2, diskTheta1, diskTheta2);
a = standardIntersection(x, yl, yr, diskRadius, t1, t2);

[x, yl, yr, t1, t2] = standardCoords(x2, x3, diskTheta1, diskTheta2);
a = a + standardIntersection(x, yl, yr, diskRadius, t1, t2);

[x, yl, yr, t1, t2]  = standardCoords(x3, x1, diskTheta1, diskTheta2);
a = a + standardIntersection(x, yl, yr, diskRadius, t1, t2);

% no orientation for final result
a = abs(a);

end

function [x, yl, yr, t1, t2] = standardCoords(x1, x2, theta1, theta2)
   % x1, x2 are coordinates, center of disk is at (0,0)
   % sign of x determines orientation.
   if norm(x1 - x2) < eps || isequal(x1, [0,0]) || isequal(x2, [0,0]) || theta1 == theta2
      yl = 0; yr = 0;
      x = 0; t1 = 0; t2 = 0;
   else
      nv = (x2-x1) / norm(x2-x1);
      yl = x1 * nv'; yr = x2 * nv';
      x = x1(1) * nv(2) - x1(2) * nv(1);

      % shift angle to open
      if x >= 0
         th = atan2(-nv(1), nv(2));
      else
         th = - atan2(-nv(1), nv(2));
      end
      t1 = theta1 - th;
      t2 = theta2 - th;
      if t1 < -pi
         t1 = t1 + 2 * pi;
         t2 = t2 + 2 * pi;
      end
      if t1 > pi
         t1 = t1 - 2 * pi;
         t2 = t2 - 2 * pi;
      end

   end
end


function area = standardIntersection(x, yl, yr, r, theta1, theta2)
   if theta1 > theta2 || yl > yr || theta2 - theta1 > 2 * pi + eps
      error('internal error')
   end

   area = 0;
   s = sign(x);
   x = abs(x);
   
   if (x <= eps)
      return
   end
   


   % calculate effective xl and xr
   % note: theta1 is in [-pi, pi] and theta2 is in [-pi, 3 pi] always theta2 > theta1
   thetal = atan2(yl, x); % thetas will be between -pi/2 and +pi/2
   thetar = atan2(yr, x);

   y1 = tan(theta1) * x;
   y2 = tan(theta2) * x;

   
   if theta1 <= thetal %  -> - pi < theta1 < pi -> theta1 < theta2 < 2 pi + theta1
      if theta2 <= thetal % both below thetal -> area = 0
         return
      else % theta2 > thetal 
         if theta2 <= thetar % theta2 in segment
            yleff = yl;
            yreff = y2; %m cut off by segment
         else % theta2 above segment
            yleff = yl;
            yreff = yr;
         end
      end
   else % theta1 > thetal
      if theta1 < thetar   % theta1 in segment
         if theta2 <= thetar % theta1 and theta2 in segment
            yleff = y1;
            yreff = y2;
         else % possibly 2 segments
            theta2e = theta2 - 2 * pi;
            if theta2e > thetal
               area = standardIntersection(s * x, yl, yr, r, theta1, thetar);
               area = area + standardIntersection(s * x, yl, yr, r, thetal, theta2e);
               return
            else
               yleff = y1;
               yreff = yr;
            end
         end
      else % theta1 above segment
         theta2e = theta2 - 2 * pi;
         if  theta2e < thetal % no intersection
            return 
         else % theta2e > thetal
            if theta2e < thetar
               yleff = yl;
               yreff = y2;
            else
               yleff = yl;
               yreff = yr;
            end
         end
      end
   end
            
  
   if (x >= r) % circle segment

      thetal = atan2(yleff, x); % thetas will be between -pi/2 and +pi/2
      thetar = atan2(yreff, x);
      area = s * r^2 * (thetar - thetal)/2;

   else % three segments, circle, triangle, cricle
      yint = sqrt(r^2 - x^2);
      area = 0;

      %left circular segment
      if yleff < - yint
         thetal = atan2(yleff, x);
         if yreff < -yint
            thetar = atan2(yreff, x);
         else
            thetar = atan2(-yint, x);
         end
         area = area + r^2 * (thetar - thetal)/2;
      end

      %middle triangle
      if yleff < -yint
         y1 = -yint;
      elseif yleff > yint
         y1 = yint;
      else
         y1 = yleff;
      end

      if yreff < -yint
         y2 = -yint;
      elseif yreff > yint
         y2 = yint;
      else
         y2 = yreff;
      end
      area = area + (y2-y1)* x / 2;

      %right circular segment
      if yreff > yint
         thetar = atan2(yreff, x);
         if yleff > yint
            thetal = atan2(yleff, x);
         else
            thetal = atan2(yint, x);
         end
         area = area + r^2 * (thetar - thetal)/2;
      end

      area = s * area;

   end

end

   


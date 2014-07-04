function a = overlapDiskTriangle(diskCenter, diskRadius, triangleCorner1, triangleCorner2, triangleCorner3)
% 
% a = overlapDiskTriangle(diskCenter, diskRadius, triangleCorner1, triangleCorner2, triangleCorner3)
% 
% description:
%     calculates the overlap between a disk and a triangle
%
% input:
%    diskCenter, diskRadius     center and radius of disk
%    triangleCorner*            corners of the triangle
%
% output:
%    a                         overlap


% relative coords

x1 = triangleCorner1 - diskCenter;
x2 = triangleCorner2 - diskCenter;
x3 = triangleCorner3 - diskCenter;


% overlaps for all 

[xl, xr, y] = standardCoords(x1, x2, diskRadius);
a = signedStandardIntersection(xl, xr, y, diskRadius);

[xl, xr, y] = standardCoords(x2, x3, diskRadius);
a = a + signedStandardIntersection(xl, xr, y, diskRadius);


[xl, xr, y] = standardCoords(x3, x1, diskRadius);
a = a + signedStandardIntersection(xl, xr, y, diskRadius);

% no orientation for final result
a = abs(a);

end

function [xl, xr, y] = standardCoords(x1, x2, r)
   % x1, x2 are coordinates, center of disk is at (0,0) 
   if isequal(x1,x2) || isequal(x1, [0,0]) || isequal(x2, [0,0])
      xl = 0; xr = 0;
      y = 2*r; % just bigger than r

   else      
      nx = (x2-x1) / norm(x2-x1);
      xl = x1 * nx'; xr = x2 * nx';
      y = x1(1) * nx(2) - x1(2) * nx(1);
   end
   
   if y < 0
      x0 = xr; xr = xl; xl = x0;
      y = -y; 
   end
   
end


function area = signedStandardIntersection(xl, xr, y, r)
   if xl > xr
      area = - standardIntersection(xr, xl, y, r);
   else
      area = standardIntersection(xl, xr, y, r);
   end
end


function area = standardIntersection(xl, xr, y, r)
      if (y <= eps)
         area = 0;
      elseif (y >= r) % circle segment
         thetal = atan2(xl, y);
         thetar = atan2(xr, y);
         area = r^2 * (thetar - thetal)/2;
      else % three segments, circle, triangle, cricle
         xint = sqrt(r^2 - y^2);  
         area = 0;
         
         %left circular segment
         if xl < - xint
            thetal = atan2(xl, y);
            if xr < -xint
               thetar = atan2(xr, y);
            else
               thetar = atan2(-xint, y);
            end
            area = area + r^2 * (thetar - thetal)/2;
         end
         
         %middle triangle
         if xl < -xint
            x1 = -xint;
         elseif xl > xint
            x1 = xint;
         else
            x1 = xl;
         end
         
         if xr < -xint
            x2 = -xint;
         elseif xr > xint
            x2 = xint;
         else
            x2 = xr;
         end 
         area = area + (x2-x1)* y / 2;
               
         %right circular segment
         if xr > xint
            thetar = atan2(xr, y);
            if xl > xint
               thetal = atan2(xl, y);
            else
               thetal = atan2(xint, y);
            end
            area = area + r^2 * (thetar - thetal)/2;
         end

      end
   
end

   


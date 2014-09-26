function pol = dilatePolygon(pol, bwidth)
%
% pol =dilatePolygon(pol, bwidth)
%
% description:
%    increase a polygon by adding a border of size bwidth
%
% input:
%    pol     array of polygon coordinates as row vectors
%    bwidth  width of the border to add (otr substract)
%
% ouput:
%    pol     polygon with added border
%
% note: for complicated shapes this simple algorithm can cause intersections !


x = pol(1,:);
y = pol(2,:);

dx = diff(x);
dy = diff(y);

dx = [dx(end), dx, dx(1)];
dy = [dy(end), dy, dy(1)];

t = 90/360 * 2 * pi;
rot = [cos(t), sin(t); -sin(t), cos(t)];

pol = zeros(2, length(x));

for i = 1:length(dx)-1  
   % find vector to push node i+1
   v1 = [dx(i), dy(i)];
   v2 = [dx(i+1), dy(i+1)]; 

   v1 = v1 / norm(v1);
   v2 = v2 / norm(v2);

   v = v1 + v2;
   v = v(:) / norm(v);
  
   vn = rot * v;

   pol(:,i) = [x(i); y(i)] + vn * bwidth;
end

end
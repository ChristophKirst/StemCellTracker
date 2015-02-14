function varargout = polygonPlot(pol, varargin)
%
% h = polygonPlot(pol, varargin)
%
% description: 
%    plots a polygon
%
% input:
%    pol    polygon
%
% output:
%    h      handle to patch object

param = parseParameter(varargin);

tri = polygonToTriangulation(pol);

h1 = patch('Faces', tri.ConnectivityList, 'Vertices', tri.Points,  'FaceColor', 'r', varargin{:}, 'EdgeColor', 'none', 'LineStyle', 'none'); 

ec = getParameter(param, 'EdgeColor', 'none');
if ~strcmp(ec, 'none')
   h = ishold;
   hold on
   
   fe = freeBoundary(tri);
   x = tri.Points(:,1); y = tri.Points(:,2);
   
   ls = getParameter(param, 'LineStyle', '-');
   lw = getParameter(param, 'LineWidth', 0.5);

   h2 = plot(x(fe)',y(fe)', [ls, ec], 'LineWidth', lw) ;

   if ~h
      hold off
   end 
end


if nargout> 0
   varargout{1} = h1;
end
if nargout> 1
   varargout{2} = h2;
end

end
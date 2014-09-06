function cl = imcolorlist(varargin)
%
% cl = imcolorlist(ncols)
% 
% descrption: 
%  returns default coloring of multi channel images
%
cl = {'r', 'g', 'b', 'w', 'm', 'y'};

if nargin > 0
   cl = padright(cl, varargin{1}, cl);
end

end
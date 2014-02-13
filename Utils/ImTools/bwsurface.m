function surface = bwsurface(bw)
%
% surface = bwsurface(bw)
%
% description: 
%    returns pixel that touch the exterior
%
% input:
%    bw    the bw labels reagion

%surface = bwdist(~bw, 'chessboard') == 1;
surface = bwdist(~bw, 'cityblock') == 1;
% faster less accurate
%surface = imfilter(double(bw), ker); % imfilter uses 0 for padding by default
%surface(surface <= conn -2) = 0;
%surface(surface > 0) = 1;

end
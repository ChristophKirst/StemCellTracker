function gr = imgrid(si)
%
% gr = imgrid(si)
%
% description:
%    returns grid with coordinates in pql format
%
% input:
%   si   size of the gird
%
% output:
%   gr   cell with grids for each coordinate, in total dim arrays

dim = length(si);

% if dim >= 2
%    si(1:2) = si([2,1]);
% end

for d = dim:-1:1
   rngs{d} = 1:si(d);
end

gr = cell(1,dim);

[gr{:}] = ndgrid(rngs{:});

% if dim >= 2
%    gr(1:2) = gr([2,1]);
%    gr{2} = flip(gr{2}, 1);
% end

end
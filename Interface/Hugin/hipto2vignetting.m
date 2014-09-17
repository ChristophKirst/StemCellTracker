function vi = hipto2vignetting(isize, pto)
%
% description:
%  create vignetting field from parsed pto file


xy = [pto(1).Vx, pto(1).Vy];

abcd = [pto(1).Va, pto(1).Vb, pto(1).Vc, pto(1).Vd];
e = pto(1).Eev;

r0 = isize/2 + xy;

rscale = 1/sqrt(sum(isize.^2))/2;

[x,y] = meshgrid(1:isize(1), 1:isize(2));

x = x - r0(1);
y = y - r0(2);

x = x * rscale;
y = y * rscale;

r2 = x.^2+y.^2;

r = ones(isize);
for i = 1:4
   vi = abcd(i) * r;
   r = r .* r2;
end

%ep = 1/2^e;
%vi = ep * vi;

end
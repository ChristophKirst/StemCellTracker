clc

t = rand(1500,25000,3);
r = rand(1500,25000,3);
dd = rand(1000,8000,3);
pos = [200,300,1];


tic

for i = 1:10
s = size(dd);
sub = {(1:s(1)) + pos(1) - 1, (1:s(2)) + pos(2) - 1,  (1:s(3)) + pos(3) - 1};
t(sub{:}) = dd;
end

toc


tic

for i = 1:10
mask = ones(size(dd));
mask = padarray(mask, pos - 1, 'pre');
mask = padarray(mask, size(r) - pos - size(dd) + 1, 'post');
mask = mask >0;

r(mask) = dd;
end
toc
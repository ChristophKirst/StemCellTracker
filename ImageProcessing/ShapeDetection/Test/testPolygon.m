

x = [1,8,3];
y = [2,7,16];

mm = poly2mask(y,x, 10,20);

mm = imdilate(mm, strel('diamond', 1)); 
ids = sub2ind(size(mm), x,y);
mm(ids) = 1;


img = zeros(10,20);
img(ids) = 1;
figure(1); clf
implot(img)


implottiling({img; mm})





%%

p1=rand(n,2);% make random polygon (self-intersecting)
pout=[sin(0:(2*pi/10) : (2*pi)); cos(0:(2*pi/10) : (2*pi))]';
pin = 0.5* p1int;

%%

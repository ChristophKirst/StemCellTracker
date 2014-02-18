img = syntheticConfocalData();



figure(52)
clf
imshow3d(img)


bw = img > 10;
figure(53)
clf
implot3d(bw)


% figure(54)
% isosurface(x,y,z,bw,0.5), 
% axis equal, title('BW')
% xlabel x, ylabel y, zlabel z
% xlim(lims), ylim(lims), zlim(lims)
% view(3), camlight, lighting gouraud

%%


3+6



%%



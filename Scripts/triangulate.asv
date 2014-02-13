
%%

idx = find(objc);
[x,y,z] = ind2sub(size(objc), idx);

tri = delaunayn([x,y,z]);

figure(50)
trisurf(tri,x,y,z);


%trisurf(tri, X(:,1), X(:,2), X(:,3))




% put single object into imaris and create surface from there ?




%% 

istack = imarisget();
size(istack)

label = istack > 0;

[v, f, nrm] = imisosurface(label);


imarisput


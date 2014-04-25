%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Imaris Interface %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test Imaris Statistics 
imaris = imarisinstance;
isurf = imarisgetobject(imaris, 'Segmentation')


%%

istat = isurf.GetStatistics

%%
%istat.mNames
istat.mFactorNames

%%

istat.mFactors

%%

rr = uint8(255 * rand(100, 100, 10));
imarissetdataset('uint8', 100, 100, 10, 1, 1);

imarissetvolume(rr)


%%

vVertices = [0,0,0;1,0,0;1,1,0];
vTriangles = [0,1,2];
vNormals = [0,0,1;0,0,1;0,0,1];
vTimeIndex = 0;
imarissetsurface(vVertices, vTriangles, vNormals);




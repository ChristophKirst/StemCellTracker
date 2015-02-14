%% Test Polygon Package

clear all
close all
clc

%% Simple Shape
clc

n = 10;
pol= 3 * [sin(0:(2*pi/n) : (2*pi)); cos(0:(2*pi/n) : (2*pi))];

%% Dilate / buffering
polb = polygonBuffer(pol, 2);
polb = polb{1};

figure(1); clf
fill(polb(1,:)',polb(2,:)','b'); hold on
fill(pol(1,:)',pol(2,:)','r')

%% Union

polb = polygonUnion({[0,0; 0,1; 1,0]', [0,0; 0,1; -1,0]'});
polb = polb{1}

figure(1); clf
fill(polb(1,:)',polb(2,:)','b'); hold on

polygonArea(polb)


%% Polygon to Triangulation
clc
tri = polygonToTriangulation(pol)

figure(2); clf; hold on
fill(pol(1,:)',pol(2,:)','r')
triplot(tri);

%% 
clc

hol =  0.3 * pol + 1;
hol2 =  0.3 * pol -1;
pol2 = 0.1 * pol + 1;


tri = polygonToTriangulation({pol, hol, hol2, pol2}, 'all', false);

figure(3); clf; hold on
fill(pol(1,:)',pol(2,:)','r');
fill(hol(1,:)',hol(2,:)','w'); fill(hol2(1,:)', hol2(2,:)', 'w');
fill(pol2(1,:)', pol2(2,:)', 'r'); 
triplot(tri);

figure(4); clf;
%patch('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'FaceColor', 'r', 'EdgeColor', 'r' )
polygonPlot({pol, hol, hol2, pol2}, 'FaceColor', 'b', 'EdgeColor', 'r', 'LineWidth', 2 )


%% Triangulation To Polygon

polt = polygonFromTriangulation(tri)

figure(5); clf;
polygonPlot(polt, 'FaceColor', 'none', 'EdgeColor', 'k')


%%


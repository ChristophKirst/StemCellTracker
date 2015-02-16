%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Clipper Interface %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Create Shape

n=7;
n2 = 13;
scale=100;

%pgon=rand(n,2);% make random polygon (self-intersecting)
pgon= 3* [sin(0:(2*pi/n) : (2*pi)); cos(0:(2*pi/n) : (2*pi))]';
phole = 2 * [sin(0:(2*pi/n2) : (2*pi)); cos(0:(2*pi/n2) : (2*pi))]';
phole = phole(:, end:-1:1);

% recast into int64 in order to input to clipper
pgonInt =int64(pgon*scale);
pholeInt = int64(phole*scale);

figure(1); clf; hold on
fill(pgonInt(:,1),pgonInt(:,2),'r'); 
fill(pholeInt(:,1),pholeInt(:,2),'w')


%% Buffer Polygon
pbuffer = mexPolygonBuffer({pgonInt', pholeInt'}, scale * 0.5, 2);
length(pbuffer)
pbgon = pbuffer{1};
pbgon = double(pbgon)' / scale;

if length(pbuffer) > 1
   pbhole = pbuffer{2};
   pbhole = double(pbhole)' / scale;
end

figure(2); clf
fill(pbgon(:,1),pbgon(:,2),'r'); hold on% plot first polygon
fill(pgon(:,1),pgon(:,2),'b')% plot first polygon

if length(pbuffer) > 1
   fill(phole(:,1), phole(:,2),'g')% plot first polygon
   fill(pbhole(:,1),pbhole(:,2),'w')% plot first polygon
end



%% Execute Operations on Polygon
clc

pol = {pgonInt', pholeInt'};
pol2 = {int64([-300,-200; -100, -100; -100, 400]')};

figure(1); clf; hold on
polygonPlot(cellfunc(@double, pol));
polygonPlot(cellfunc(@double, pol2), 'FaceColor', 'b')


%%

pe = mexPolygonExecute(1, 0, pol, pol2);
length(pe)

figure(2); clf
polygonPlot(pe)


%% Plotting

%patch('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'FaceColor', 'r', 'EdgeColor', 'r' )



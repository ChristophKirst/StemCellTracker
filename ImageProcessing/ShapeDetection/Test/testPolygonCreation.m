%%
load ./ImageProcessing/ShapeDetection/Test/points.mat
figure(1); clf
[v, s] = detectAlphaVolume(pp, 100, true);


%%
T = unique(s.tri, 'rows');
t = triangulation(T, pp);

figure(7); clf;
triplot(t)

%%
% construct triangle connectivity graph
n = size(pp,1)
mat = sparse(n,n);

ea = t.edgeAttachments(t.edges)
ea = ea(cellfun(@(x) length(x) > 1, ea))

% in 2d only 2 triangles can share an edge
A = graphEdgesToAdjacencyMatrix(cell2mat(ea));
cc = graphAdjacencyMatrixToConnectedComponents(A)


% convert free boundaries to borders


%%
nc =length(cc);
col = colorcube(nc+6);


figure(7); clf; hold on
for i = 1:nc
   trij = triangulation(T(cc{i}, :), pp);
   triplot(trij, 'color', col(i,:))
end

%%


pol = [sin(0:(2*pi)/5:(2*pi)); cos(0:(2*pi)/5:(2*pi))];
X = pol(1,:)'; Y= pol(2,:)';
figure(7); clf;
plot(X,Y, '+')

%%
clc
polb = bufferPolygon(X,Y,2, 'out')



%%









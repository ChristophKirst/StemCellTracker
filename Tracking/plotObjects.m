function plotObjects(data)
%
% plotObjects(data)
%
% description:
%    plots the objects in array data
% 

xyz = [data.r];
s = [data.volume];
s = s / max(s) .* 200;
col = (1:size(xyz,2)) / size(xyz,2);

if data.dim ==2
   scatter(xyz(1,:), xyz(2,:), s, col)
else
   scatter3(xyz(1,:), xyz(2,:), xyz(3,:), s, col)
end
   
end


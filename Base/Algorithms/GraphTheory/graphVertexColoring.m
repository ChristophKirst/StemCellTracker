function cids = graphVertexColoring(v, e)
%
% cids = graphVertexColoring(v, e)
%
% description:
%    colors a graph (v,e) using the greedy non-optimal dsatur algorithm 
%
% input:
%    v   vertices
%    e   edges as rows [node_pre, node_post; ...]
%
% output:
%    cids  coloring ids

v = v(:);

if isempty(e)
   cids = ones(size(v,1),1);
   return
end


n = size(v,1);
cids = zeros(n,1);
available_colors = 1;

% relabel vertices and eges to 1:n
[~, rel] = ismember(1:max(v), v);
e = rel(e);

% degrees
for i=n:-1:1
    v = i;
    degrees(i,1) = size([e(e(:,1)==v,2); e(e(:,2)==v,1)],1);
end

% degrees of saturation
degrees_of_saturation = zeros(n,1);

% Coloring
for i=1:n
    if i == 1
        [~, index] = max(degrees);
        v = index(1);
        cids(v) = 1;
        assigned_color_v = 1;
    else
       uncolored = find(cids==0);
       index_temp = degrees_of_saturation(uncolored)==max(degrees_of_saturation(uncolored));
       index = uncolored(index_temp);
       if(size(index,1)>1)
          [~, index1] = max(degrees(index));
          v = index(index1(1));
       else
          v = index;
       end
       
       % Assign first available color to v
       neighbors = [e(e(:,1)==v,2); e(e(:,2)==v,1)];
       for j=1:available_colors
          if size(find(cids(neighbors)==j),1)==0
             cids(v) = j;
             assigned_color_v = j;
             break;
          end
       end
       if cids(v) == 0
          available_colors = available_colors + 1;
          cids(v) = available_colors;
          assigned_color_v = available_colors;
       end
    end

    % Update Degrees of saturation
    neighbors_v = [e(e(:,1)==v,2); e(e(:,2)==v,1)];
    for j=1:size(neighbors_v,1)
       u = neighbors_v(j);
       neighbors_u = [e(e(:,1)==u,2); e(e(:,2)==u,1)];
       if size(find(cids(neighbors_u)==assigned_color_v),1) == 1
           degrees_of_saturation(u,1) = degrees_of_saturation(u,1) + 1;
       end
    end

end
function a = edges2AdjacencyMatrix(e, varargin)
%
% a = connectedAlignments(e, n)
%
% description:
%    converts edges e to adjacenty matrix a of a graph 
%
% input:
%    e     edges as array of pairs [pair1; pair2;...]   absolute shifts 
%    n     (optional) number of nodes (max(e(:)))
%
% output:
%    a     adjacency matrix 

if nargin > 1
   n = varargin{1};
else
   n = max(e(:));
end

a = sparse(e(:,1), e(:,2), 1, n, n);

end
 
   
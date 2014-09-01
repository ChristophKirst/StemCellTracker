function c = adjacencyMatrix2ConnectedComponents(a)
%
% c = adjacencyMatrix2ConnectedComponents(a)
% 
% description:
%    returns indices of nodes as cell array of the conencted components of a graph
%    represented by the symmetric adjacency matrix a
%
% input:
%    a   adjacency matrix
%
% output:
%    c   connected components as cell array of indices arrays
%
% See also: dmperm


% Check size of adjacency matrix
[n, m] = size(a);
if n ~= m
   error ('adjacencyMatrix2ConnectedComponents: Adjacency matrix must be square!')
end;

% make shure matrix is symmetric
a = a + a';
a = a > 0;

% Dulmage-Mendelsohn permutation on A
if ~all(diag(a))
   [~, p, ~, r] = dmperm(a | speye(size(a)));
else
   [~, p, ~, r] = dmperm(a);
end

% cell array of indices

c = cell(1,length(r)-1);
for i = 1:length(c)
   c{i} = p(r(i):r(i+1)-1);
end

end
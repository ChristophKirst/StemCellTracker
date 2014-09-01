function comp = connectedAlignments(a, varargin)
%
% comp = connectedAlignments(a, param)
%
% description:
%    uses the qulaity measure in the AlignmentPairs in Alignment a to determine
%    connected components and returns these as a cell array of Alignment classes
%
% input:
%    a         Alignment class
%    param     parameter struct with entries
%              .threshold.quality     quality threshold (-Inf)
%
% output:
%    comp      conencted components as cell araay of Alignment classes
%
% See also: Alignment

if ~isa(alignm, 'Alignment')
   error('connectedAlignments: expects Alignment class as input');
end

param = parseParameter(varargin{:});
thq = getParameter(param, 'threshold.quality', -Inf);

if thq == -Inf
   comp = {a};
   return
end

% construct adjacency matrix and find connected components

e = [];
for p = 1:alignm.npairs

   if a.pairs(p).quality > thq
      e = [e; [a.pairs(p).from, a.pairs(p).to]];
   end
end

adj = edges2AdjacencyMatrix(e);
c = adjacencyMatrix2ConnectedComponents(adj);


% construc Alignment classes

comp = cell(1, length(c));





end
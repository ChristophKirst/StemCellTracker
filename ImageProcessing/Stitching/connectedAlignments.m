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
%    comp      connected components as array of Alignment classes
%
% See also: Alignment

if ~isa(a, 'Alignment') %%&& ~isa(a,'ImageSourceAligned')
   error('connectedAlignments: expects Alignment class as input');
end

param = parseParameter(varargin{:});
thq = getParameter(param, 'threshold.quality', -Inf);

if thq == -Inf
   comp = a;
   return
end

% construct adjacency matrix and find connected components

e = [];
pairs = a.pairs;

for p = 1:a.nPairs
   if pairs(p).quality > thq
      e = [e; [pairs(p).from, pairs(p).to]]; %#ok<AGROW>
   end
end

adj = graphEdgesToAdjacencyMatrix(e, a.nNodes);
c = graphAdjacencyMatrixToConnectedComponents(adj);

% construct Alignment / ImageSourceAligned classes

nodes = a.nodes;

for i = length(c):-1:1
   if isa(a, 'ImageSourceAligned')
      as = ImageSourceAligned(a);
   else
      as = Alignment(a);
   end
   
   as.anodes = nodes(c{i});
   as.reducePairs();
   as.removeLowQualityPairs(thq);

   comp(i) = as;
end

end
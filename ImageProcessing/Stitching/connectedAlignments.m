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

%dataSize = a.dataSize;
%overlap = getParameter(param, 'overlap', 0);
%overlap = padright(overlap, length(dataSize), overlap);

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

o = a.position;
nodes = a.nodes;

for i = length(c):-1:1 
   as = Alignment(a);

   as.anodes = nodes(c{i});
   as.reducePairs();
   as.removeLowQualityPairs(thq);
   
   % estimate new origin:
   %onode = as.anodes(1);
   %sub = a.source.cellIndexToSubIndex(onode);
   %as.aorigin  = (dataSize - overlap) .* (sub -1) + o;
   
   comp(i) = as;
end

end
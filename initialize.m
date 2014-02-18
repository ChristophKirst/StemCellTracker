function [] = initialize(comp)
%
% [] = initialize()
%
% description:
%     sets all necessary paths and starts up interfaces
%
% See also: setPath

% paths
setPath();
     
% figure output on first screen
set(0, 'DefaultFigurePosition',  [1   705   560   420]);

if nargin < 1
   comp = 0;
end

if comp
   cd Segmentation
   mex segmentByPropagationMEX.cpp
   cd ..
end
   


% imagej interface
% if nargin < 1
%    with_ij = 0;
% end
% if with_ij
%    ijstart
% end


end


function [] = initialize(comp)
%
% [] = initialize()
%
% description:
%     sets all necessary paths and parameter
%
% See also: setPath

% paths
setPath();
     
% figure output on first screen
set(0, 'DefaultFigurePosition',  [1   705   560   420]);


% compile code
if nargin < 1
   comp = 0;
end

if comp
   compileSegmentation();
   ijcompile();
end

end


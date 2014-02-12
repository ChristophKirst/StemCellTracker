function [] = initialize(with_ij)
%
% [] = initialize()
%
% description:
%     sets all necessary paths etc
%

% paths
addpath('./Classes', ...
        './Filtering', ...
        './Segmentation', ...
        './Tracking',...
        './Interface',...
        './Interface/ImageFormats',...
        './Utils',...
        './Utils/ImTools',...
        './Utils/ImageJ',...
        './Utils/Imaris');
 
 addpath('./Test/Scripts');
        

     
% figure output on first screen
set(0, 'DefaultFigurePosition',  [1   705   560   420]);


% imagej interface
if nargin < 1
   with_ij = 0;
end

if with_ij
   ijstart
end

end


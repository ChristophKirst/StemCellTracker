%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ImageSourceFiles - test different folders data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image data is read from disk each time and not cached
% useful for large tilings that do not fit into memory

clear all
clear classes
close all
clc
initialize
bfinitialize

verbose = true;


%% 
clc

% infer the tag expression o he files from the file folder automatically
texp = tagExpression('./Test/Images/hESCells_Cytoo/*.tif', 'tagnames', {'S', 'C'})
%texp='140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<pos,6>t00000001z001c<ch,2>.tif'

%%
is = ImageSourceFiles(texp);
is.printInfo



%%
% infer the tag expression o he files from the file folder automatically
texp = tagExpression('./Test/Images/hESCells_GCamp_Vignetting/*.tif', 'tagnames', {'S', 'C'})
%texp='140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<pos,6>t00000001z001c<ch,2>.tif'

%%
is = ImageSourceFiles(texp);
is.printInfo

is

%%

is.resetRange;
is.setReshape('S', 'UV', [3,3]);
is.setRawCellDataCaching(true);
is.setCellFormat('Uv');

is.printInfo

%%
is.range
is.cellSize
is.cellFormat


%%
is.plottiling


%% restric range to some sub set

is.resetRange;
is.cellIndex

%%
is.cellSize

%%
is.setRange('U', 1:2, 'V', 2:3);
is.cellFormat
is.cellSize
is.rawCellIndex

%%

is.rangeFromVarargin


%%
is.cellIndex
is.cellSize

%%

is.plottiling


%%
clc
is.resetRange;
is.cellIndex('U', 1:2, 'v', 1:2)


%%

% % check formatting
% cd = is.cell(1:18)
% figure(1); clf
% implottiling(cd, 'tiling', [8,2])

 
%%

% make preview
preview = is.cellPreview('overlap', 110, 'scale', 0.1, 'lines', true);
figure(2); clf
implot(preview)


%%

figure(3)
%is.plottiling % this crashes for large tilings !















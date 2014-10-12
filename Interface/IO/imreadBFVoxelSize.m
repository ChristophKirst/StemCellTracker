function vsize = imreadBFVoxelSize(name, varargin)
%
% vsize = imreadBFVoxelSize(name, varargin)
%
% decription:
%    tries to infer voxel sizes from bio-format image
%
% input:
%    name    filename or reader
%    param   parameter struct with entries
%            .S / .s  series ids
%             
% output:
%    emwl    emission wave length
%    exwl    excitation wave length
%

param = parseParameter(varargin);

ireader = imreadBFReader(name);

% series specifed

[sids, nSeries] = imfrmtParseRangeToIndex(ireader.getSeriesCount(), 'S', param);
  
metadataStore = ireader.getMetadataStore();

vsize = cell(nSeries, 1);

sc = 1;
for s = sids

   % Voxel size X
   voxelX = metadataStore.getPixelsPhysicalSizeX(s-1);
   if isempty(voxelX)
      voxelX = 0;
   else
      voxelX = voxelX.getValue();
   end
   
   % Voxel size Y
   voxelY = metadataStore.getPixelsPhysicalSizeY(s-1);
   if isempty(voxelY)
      voxelY = 0;
   else
      voxelY = voxelY.getValue();
   end
   
   % Workaround for 2-D ND2 files
   if voxelY == 0 && voxelX > 0
      voxelY = voxelX;
   end
   
   % Voxel size Z
   voxelZ = metadataStore.getPixelsPhysicalSizeZ(s-1);
   if isempty(voxelZ)
      voxelZ = 0;
   else
      voxelZ = voxelZ.getValue();
   end
   voxels = [voxelX, voxelY, voxelZ]';
   
   % set data
   vsize{sc} = voxels;
   sc = sc + 1;
end
   
    
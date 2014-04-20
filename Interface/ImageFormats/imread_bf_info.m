function meta = imread_bf_info(fname)
%
% meta = imread_bf_info(name, param) 
%
% description:
%    import bio-format meta data from file using the loci tools
%
% input:
%    fname   (optional) filename
%
% output:
%    meta    meta data indicating time specifications
%
% note: the loci_tools.jar must be present in the javaclasspath
%       try bfinitialize or ijintialize to set this path
%
% See also: imread_bf, bfinitialize, imread, ijimage2mat

if nargin == 0
  [file, path] = uigetfile(bffileextensions, 'Choose a file to open');
  fname = [path file];
end
if ~exist(fname, 'file')
   error('imread_bf: file does not exists: %s', fname)
end


% initialize logging
loci.common.DebugTools.enableLogging('INFO');

% Get the channel filler
r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
%if stitchFiles
%    r = loci.formats.FileStitcher(r);
%end
r.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
r.setId(fname);

numSeries = r.getSeriesCount();

for s = numSeries:-1:1

   %%% read a single series 
   meta(s).series = s;

   r.setSeries(s - 1);
   meta(s).images = r.getImageCount();

   % dimensions
   meta(s).sizeX    = r.getSizeX();
   meta(s).sizeY    = r.getSizeY();
   meta(s).sizeZ    = r.getSizeZ();
   meta(s).channels = r.getSizeC();
   meta(s).times    = r.getSizeT();

   %meta.width    = r.getSizeX();
   %meta.height   = r.getSizeY();
   %meta.depth    = r.getSizeZ();
   %meta.frames   = r.getSizeT();
   %meta.channels = r.getSizeC();

   metadataList = r.getMetadata();
   %m=r.getMetadataStore();
   subject = metadataList.get('parameter scale');

   if ~isempty(subject)% if possible pixelsizes are added (only ics files)
      voxelsizes=str2num(subject); %#ok<ST2NM>

      meta(s).voxelsizes=voxelsizes;
   
      if voxelsizes>1
         meta(s).PixelSizeX = voxelsizes(2);
         meta(s).PixelSizeY = voxelsizes(3);
         meta(s).PixelSizeZ = voxelsizes(4);
         meta(s).PixelSizeT = voxelsizes(5);
      end
   end

   meta(s).raw=metadataList;%metadata as java hashtable

   %readout of java hashtable
   meta(s).rawNames = meta(s).raw.keySet.toArray;
   meta(s).rawValues= meta(s).raw.values.toArray;
   
   meta(s).parameterNames = cell(meta(s).raw.keySet.toArray);
   meta(s).parameterValues= cell(meta(s).raw.values.toArray);
   
   meta(s).parameters = table(meta(s).parameterNames, meta(s).parameterValues, 'VariableNames', {'Parameter', 'Value'});
end

end
            
   
        
   
   
   
   
   
            
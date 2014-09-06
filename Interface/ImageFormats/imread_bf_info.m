function iinfo = imread_bf_info(fname, varargin)
%
% iinfo = imread_bf_info(name, param) 
%
% description:
%    import bio-format meta data from file using the loci tools
%
% input:
%    fname   (optional) filename, if not given open file browser
%    param   (optional) parameter struct with entries
%            .metadata  include meta data in iinfo (true)
%            .table     transform metadata to a table (true)
%            .struct    return sturct instead of ImageInfo class (false)
%
% output:
%    iinfo   struct/ImageInfo class containing the info of the image, for multiple series an array of this class is returned
%
% note: 
%    the loci_tools.jar must be present in the javaclasspath
%    use bfinitialize to initialize the loci tools.
%
% See also: ImageInfo, imread_bf, bfinitialize, imread

if nargin == 0
  [file, path] = uigetfile(bffileextensions, 'Choose a file to open');
  fname = [path file];
end
if ~exist(fname, 'file')
   error('imread_bf_info: file does not exists: %s', fname)
end

param = parseParameter(varargin{:});
meta  = getParameter(param, 'metadata', true);
tab   = getParameter(param, 'table',    true);

% initialize logging
% loci.common.DebugTools.enableLogging('INFO');

% Get the channel filler
r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
%if stitchFiles
%    r = loci.formats.FileStitcher(r);
%end
r.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
r.setId(fname);

numSeries = r.getSeriesCount();

if ~getParameter(param, 'struct', false)
   iinfo(numSeries) = ImageInfo();
end

for s = numSeries:-1:1

   %%% read a single series 
   iinfo(s).iseries = s;

   r.setSeries(s - 1);
   iinfo(s).inimages = r.getImageCount();

   % dimensions
   iinfo(s).isizePQLCT = [r.getSizeY(), r.getSizeX(), r.getSizeZ(), r.getSizeC(), r.getSizeT()];
   
   isize = iinfo(s).isizePQLCT;
   iinfo(s).isize = isize(isize ~= 1); % size after squeezing
   
   iformat = 'pqlct';
   iinfo(s).iformat = iformat(isize ~=1); % format after squeeze
   
   iinfo(s).iclass = 'double'; % bf read always returns double !

   % meta data
   if meta
      iinfo(s).imetadata  = r.getMetadata();   %metadata as java hashtable % %m=r.getMetadataStore()
    
      iinfo(s) = determineScale(iinfo(s));
      iinfo(s) = determineColor(iinfo(s));
      
      if tab
         % readout of java hashtable
         %rawNames = iinfo(s).metadata.keySet.toArray;
         %rawValues= iinfo(s).metadata.values.toArray;
   
         parameterNames  = cell(iinfo(s).imetadata.keySet.toArray);
         parameterValues = cell(iinfo(s).imetadata.values.toArray);
   
         % make a table out of it
         iinfo(s).imetadata = table(parameterNames, parameterValues, 'VariableNames', {'Parameter', 'Value'});
      end   
   else
      iinfo(s).imetadata = [];
   end
 
end

end
            
   
        
   
% helper

%infer spatial scales from meta data, sets: scaleP, scaleQ, scaleL, times
function iinfo = determineScale(iinfo)
   % extract info from meta data:

   % TIF XYResolution is number of Pixel per ResolutionUnit
   xr = iinfo.imetadata.get('XResolution');
   if ~isempty(xr)
      yr = iinfo.imetadata.get('YResolution');

      iinfo.iscale = [yr, xr];
      %iinfo.iscaleFormat = 'pq';
      iinfo.iunit =  iinfo.imetadata.get('ResolutionUnit');

      return
   end


   %TODO:  %ZVI  %LIF
   % ICS todo: check
   voxelsizes = iinfo.imetadata.get('parameter scale');
   if ~isempty(voxelsizes)% if possible pixelsizes are added (only ics files)
      voxelsizes=str2num(voxelsizes); %#ok<ST2NM>

      if voxelsizes>1
         iinfo.iscale = voxelsizes(2:5);
         %iinfo.iscaleFormat = 'pqlt';
      end

      return
   end
  
end




%infer colors
function iinfo = determineColor(iinfo)

   % TIF PhotometricInterpretation gives color interpretation
   c = iinfo.imetadata.get('PhotometricInterpretation');
   if ~isempty(c)
      switch c
         case 'BlackIsZero'
            c = {'gray'};
         case 'RGB'
            c = {'r', 'g', 'b'};
         otherwise
            c = repmat({'gray'}, 1, iinfo.sizeC);
      end

      iinfo.icolor = c;

      return
   end


   %TODO:
   %ZVI
   %LIF
 
end
 

   
   
            
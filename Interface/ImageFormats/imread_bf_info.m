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
%            .struct    return struct instead of ImageInfo class (false)
%            .array     return image info for each series separately (false, i.e. assuming it is the same for all series, e.g. for a tiled image)
%            .squeeze   return info assuming data is squezed by imread_bf (true)
%            sub image specs as passed imread_bf to determine the returned size
%            .series      ids for seris (if array data is cell array) ([] = all)
%            .time /.t    ids of time frames to import [tid1, tid2, ...], {tmin, tmax} ([] = all)
%            .channel /.c ids of channels to import [cid1, cid2, ...] {cmin, cmax} ([] = all)
%            .z /.l       pixel ids in z / l direction {zmin, zmax} ([] = all)
%            .x /.p       pixel ids in x / p direction {xmin, xmax} ([] = all)
%            .y /.q       pixel ids in y / q direction {ymin, ymax} ([] = all)
%
% output:
%    iinfo   struct/ImageInfo class containing the info of the image, for multiple series an array of this class is returned
%
% note: 
%    the loci_tools.jar must be present in the javaclasspath
%    use bfinitialize to initialize the loci tools.
%
% See also: imread_bf, ImageInfo, bfinitialize, imread

if nargin == 0
  [file, path] = uigetfile(bffileextensions, 'Choose a file to open');
  fname = [path file];
end
if ~isfile(fname)
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

sq = getParameter(param, 'squeeze', true);
ar = getParameter(param, 'array',   false);



seriesid = getParameter(param, 'series', []);
if isempty(seriesid)
   seriesid = 1:r.getSeriesCount();
end
numSeries = length(seriesid);

if ar
   
   if ~getParameter(param, 'struct', false)
      iinfo(numSeries) = ImageInfo();
   end

   for s = seriesid
      %%% read a single series
      iinfo(s).iseries = s;
      
      r.setSeries(s - 1);
      iinfo(s).inimages = r.getImageCount();
      
      % dimensions
      isizePQLCT = [0,0,0,0,0];
      
      isizePQLCT(5) = getSize(param, 'time'   , 't', r.getSizeT());
      isizePQLCT(4) = getSize(param, 'channel', 'c', r.getSizeC());
      isizePQLCT(3) = getSize(param, 'z',       'l', r.getSizeZ());
      isizePQLCT(2) = getSize(param, 'y',       'q', r.getSizeY());
      isizePQLCT(1) = getSize(param, 'x',       'p', r.getSizeX());

      iinfo(s).idatasizePQLCT = isizePQLCT;
      
      isize = iinfo(s).idatasizePQLCT;
      iformat = 'pqlct';
   
      if sq
         iinfo(s).idatasize   = isize(isize ~= 1); % size after squeezing
         iinfo(s).idataformat = iformat(isize ~=1); % format after squeeze
      else
         iinfo(s).idatasize = isize;
         iinfo(s).idataformat = iformat;
      end

      iinfo(s).idataclass = 'double'; % bf read always returns double !
      
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
      
      iinfo(s).irawformat = iinfo(s).idataformat;
      iinfo(s).irawsize   = iinfo(s).idatasize;
      iinfo(s).icellsize = numSeries;
      iinfo(s).icellformat ='u';    
   end

else
   
   if ~getParameter(param, 'struct', false)
      iinfo = ImageInfo();
   end

   %%% single series 
   iinfo.iseries = 1:numSeries;

   r.setSeries(0);
   iinfo.inimages = r.getImageCount();

   % dimensions
   % dimensions
   isizePQLCT = [0,0,0,0,0];
   
   isizePQLCT(5) = getSize(param, 'time'   , 't', r.getSizeT());
   isizePQLCT(4) = getSize(param, 'channel', 'c', r.getSizeC());
   isizePQLCT(3) = getSize(param, 'z',       'l', r.getSizeZ());
   isizePQLCT(2) = getSize(param, 'y',       'q', r.getSizeY());
   isizePQLCT(1) = getSize(param, 'x',       'p', r.getSizeX());
   
   iinfo.idatasizePQLCT = isizePQLCT;
   
   isize = iinfo.idatasizePQLCT;
   iformat = 'pqlct';
   
   if sq
      iinfo.idatasize = isize(isize ~= 1); % size after squeezing
      iinfo.idataformat = iformat(isize ~=1); % format after squeeze
   else      
      iinfo.idatasize = isize;
      iinfo.idataformat = iformat; 
   end
   
   iinfo.idataclass = 'double'; % bf read always returns double !

   % meta data
   if meta
      iinfo.imetadata  = r.getMetadata();   %metadata as java hashtable % %m=r.getMetadataStore()
    
      iinfo = determineScale(iinfo);
      iinfo = determineColor(iinfo);
      
      if tab
         % readout of java hashtable
         %rawNames = iinfo(s).metadata.keySet.toArray;
         %rawValues= iinfo(s).metadata.values.toArray;
   
         parameterNames  = cell(iinfo.imetadata.keySet.toArray);
         parameterValues = cell(iinfo.imetadata.values.toArray);
   
         % make a table out of it
         iinfo.imetadata = table(parameterNames, parameterValues, 'VariableNames', {'Parameter', 'Value'});
      end   
   else
      iinfo.imetadata = [];
   end
   
   iinfo.irawformat = iinfo.idataformat;
   iinfo.icellsize = numSeries;
   iinfo.icellformat ='u';
 
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

      iinfo.iscale = [xr, yr];
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
            c = imcolorlist(iinfo.isizePQLCT(4));
      end

      iinfo.icolor = c;

      return
   end

   iinfo.icolor = imcolorlist(iinfo.isizePQLCT(4));
   
   %TODO:
   %ZVI
   %LIF
 
end
 

   

function si = getSize(param, name, nameshort, isize)
      
   if isentry(param, name)
      si = getParameter(param, name, []);
   else
      si = getParameter(param, nameshort, []);
   end
   

   if isempty(si)
      si = isize;
   elseif iscell(si)
%       if length(si) ~= 2
%          error(['imread_bf_info: ', name, ' indices not a cell of {minid,maxid} but %s!'], var2char(si))
%       end
%       
%       if isempty(si{2}) || ischar(si{2}) || si{2} > isize  % char = 'end'
%          si{2} = isize;
%       end
%       if isempty(si{1}) || ischar(si{1}) || si{1} < 1
%          si{1} = 1;
%       end
%       si = si{2} - si{1};
      si = cell2mat(si);
   end
 
   if si > isize || si < 1
      error(['imread_bf_info: ', name, ' size out of range !'])
   end
end
   
            
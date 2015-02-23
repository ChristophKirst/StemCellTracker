function [iinfo, ireader] = imreadBFInfo(fname, varargin)
%
% iinfo = imreadBFInfo(name, param) 
% [iinfo, ireader] = imreadBFInfo(name, param) 
% [iinfo, ireader] = imreadBFInfo(name, iinfo, param) 
%
% description:
%    import bio-format meta data from file using the loci tools
%
% input:
%    fname   (optional) filename, ireader, if empty open file browser
%    iinfo   (optional) ImageInfo class, only the raw data info is over wrtten in this case and the data sizes are updated
%    param   (optional) parameter struct with entries
%            .metadata  include meta data in iinfo (true)
%            .table     transform metadata to a table (true)
%            .struct    return struct instead of ImageInfo class (false)
%            .array     return image info for each series separately (false, i.e. assuming it is the same for all series, e.g. tiles of a larger image)
%            .squeeze   return info assuming data is squeezed by imreadBF (true)
%
%            sub image specs as passed to imreadBF to determine the returned size:
%            .S / .s   ids for series     [sid1, sid2, ...], {sid1, sid2, ...} ([] = all)
%            .T / .t   ids of time frames [tid1, tid2, ...], {tid1, tid2, ...} ([] = all)
%            .C / .c   ids of channels    [cid1, cid2, ...]. {cid1, cid2, ...} ([] = all)
%            .Z / .z   pixel ids in z     [zid1, zid2, ...], {zid1, zid2, ...} ([] = all)
%            .X / .x   pixel ids in x     [xid1, xid2, ...], {xid1, xid2, ...} ([] = all)
%            .Y / .y   pixel ids in y     [yid1, yid2, ...], {yid1, yid2, ...} ([] = all)
%
% output:
%    iinfo   struct/ImageInfo class containing the info of the image
%            for multiple series an array of this class is returned
%    ireader (optional) channel reader java class for subsequent image access
%
% note: 
%    the loci_tools.jar must be present in the javaclasspath
%    use bfinitialize to initialize the loci tools.
%
% See also: imreadBF, ImageInfo, bfinitialize, imread

ireader = imreadBFReader(fname);

fromInfo = false;
if nargin > 2
   if isa(varargin{1}, 'ImageInfo')
      infoIn = varargin{1};
      fromInfo = true;
      varargin = varargin{2:end};
   end
end

param = parseParameter(varargin{:});

meta  = getParameter(param, 'metadata', true);
tab   = getParameter(param, 'table',    true);

sq = getParameter(param, 'squeeze', true);
ar = getParameter(param, 'array',   false);

[sids, S] = imfrmtParseRangeToIndex(ireader.getSeriesCount(), 'S', param);

if ar  
   if ~getParameter(param, 'struct', false)
      if fromInfo
         iinfo(S) = infoIn.copy();
      else
         iinfo(S) = ImageInfo();
      end
   end

   % loop over all series
   for s = sids
      
      %%% read a single series
      ireader.setSeries(s - 1);
      
      % dimensions
      isize = [0,0,0,0,0];
      isize(5) = getSize(param, 'T', 't', r.getSizeT());
      isize(4) = getSize(param, 'C', 'c', r.getSizeC());
      isize(3) = getSize(param, 'Z', 'z', r.getSizeZ());
      isize(2) = getSize(param, 'Y', 'y', r.getSizeY());
      isize(1) = getSize(param, 'X', 'x', r.getSizeX());

      % format 
      iformat = 'XYZCT';

      if sq
         iinfo(s).irawdatasize   = isize(isize ~= 1);  % size after squeezing
         iinfo(s).irawdataformat = iformat(isize ~=1); % format after squeeze
      else
         iinfo(s).irawdatasize = isize;
         iinfo(s).irawdataformat = iformat;
      end

      iinfo(s).irawdataclass  = 'double'; % bf read always returns double 
      
      iinfo(s).irawcellsize   = 1;
      iinfo(s).irawcellformat ='S';
      
      % set raw specs if ImageInfo class
      if isa(iinfo, 'ImageInfo')  
         iinfo(s).initializeCellDataFormatsAndSizesFromRaw;
         iinfo(s).initializeDataAndCellSizeFromRaw;
         iinfo(s).initializeDataClassFromRaw;
      end

      % meta data
      if meta
         
         % load meta data
         iinfo(s).imetadata  = r.getMetadata();   %metadata as java hashtable % 
         %m=r.getMetadataStore()
         %seriesMetadata = r.getSeriesMetadata();
         %loci.formats.MetadataTools.merge(globalMetadata, seriesMetadata, 'Global ');
         
         % extract scales
         iinfo(s) = determineScale(ireader, iinfo(s), s);
         iinfo(s) = determineColor(ireader, iinfo(s), s);
         iinfo(s) = determineChannelname(ireader, iinfo(s), s);
         
         % convert to table
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

else % read info from first spcified series only
   
   if ~getParameter(param, 'struct', false)
      if fromInfo
         iinfo = infoIn.copy();
      else
         iinfo = ImageInfo();
      end
   end

   % series 
   ireader.setSeries(sids(1)-1);

   % dimensions
   isize = [0,0,0,0,0];
   isize(5) = getSize(param, 'T', 't', ireader.getSizeT());
   isize(4) = getSize(param, 'C', 'c', ireader.getSizeC());
   isize(3) = getSize(param, 'Z', 'z', ireader.getSizeZ());
   isize(2) = getSize(param, 'Y', 'y', ireader.getSizeY());
   isize(1) = getSize(param, 'X', 'x', ireader.getSizeX());
   
   % format
   iformat = 'XYZCT';
   if sq
      iinfo.irawdatasize   = isize(isize ~= 1);  % size after squeezing
      iinfo.irawdataformat = iformat(isize ~=1); % format after squeeze
   else
      iinfo.irawdatasize = isize;
      iinfo.irawdataformat = iformat;
   end
   
   iinfo.irawdataclass  = 'double'; % bf read always returns double
   
   iinfo.irawcellsize   = S;
   iinfo.irawcellformat ='S';
 
   % set raw specs if ImageInfo class
   if isa(iinfo, 'ImageInfo')
      iinfo.initializeCellDataSizeAndFormatFromRaw;
   end

   
   % meta data
   if meta
      
      % load meta data
      iinfo.imetadata  = ireader.getMetadata();   %metadata as java hashtable % %m=r.getMetadataStore()
    
      % detremined properties
      iinfo = determineScale(ireader, iinfo, 1);
      iinfo = determineColor(ireader, iinfo, 1);
      iinfo = determineChannelname(ireader, iinfo, 1);
      
      % meta data to table
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
end

% close reader if not returned
if nargout <= 1
   ireader.close();
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper

% parse ids
function si = getSize(param, name, nameshort, isize)
      
   if isentry(param, name)
      si = getParameter(param, name, []);
   else
      si = getParameter(param, nameshort, []);
   end

   if isempty(si)
      si = isize;
   elseif iscell(si)
      si = cell2mat(si);
   end
 
   if si > isize || si < 1
      error(['imreadBFInfo: ', name, ' size out of range !'])
   end
end



%infer spatial scales from meta data, sets: scaleP, scaleQ, scaleL, times
function iinfo = determineScale(ireader, iinfo, s)
   vs = imreadBFVoxelSize(ireader, 'S', s);
   vs = vs{1};
   
   iinfo.setScale('X', vs(1));
   iinfo.setScale('Y', vs(2));
   iinfo.setScale('Z', vs(3));
end
   


%infer colors
function iinfo = determineColor(~, iinfo, ~)
   % TIF PhotometricInterpretation gives color interpretation
   c = iinfo.imetadata.get('PhotometricInterpretation');
   if ~isempty(c)
      switch c
         case 'BlackIsZero'
            c = {'gray'};
         case 'RGB'
            c = {'r', 'g', 'b'};
         otherwise
            c = imcolorlist(iinfo.datasizeC);
      end

      iinfo.icolor = c;

      return
   end

   iinfo.icolor = imcolorlist(iinfo.dataSizeC);
   
   %TODO:
   %ZVI
   %LIF
end
  


%infer channel name
function iinfo = determineChannelname(ireader, iinfo, s)
   cn = imreadBFChannelNames(ireader, 'S', s);
   if all(cellfun(@length, cn))
      iinfo.setChannelName(cn(:,1)');
   end
end

















% 
% % Extract the color for all channels
%     function colors = getChannelColors(ireader)
%         
%         % Get the metadata store
%         metadataStore = ireader.getMetadataStore();
%         
%         nColors = ireader.getSizeC();
%         colors = zeros(nColors, 4);
%         
%         for ch = 1 : nColors
%             
%             try
%                 color = metadataStore.getChannelColor(currentSeries, ch - 1);
%             catch
%                 color = [];
%             end
%             if isempty(color)
%                 colors = [255 255 255 255];
%             else
%                 colors(ch, 1 : 4) = [ ...
%                     color.getRed(), color.getGreen(), ...
%                     color.getBlue(), color.getAlpha()];
%             end
%             
%         end
%         
%         % Return colors in the 0 .. 1 range
%         if max(colors(:)) <= 255
%             colors = colors ./ 255;
%         end
%         
%     end
% 


% %%TODO: use imreadBFVoxelSize and imreadBFTimeStamps
% 
%    % extract info from meta data:
% 
%    % TIF XYResolution is number of Pixel per ResolutionUnit
%    xr = iinfo.imetadata.get('XResolution');
%    if ~isempty(xr)
%       yr = iinfo.imetadata.get('YResolution');
% 
%       iinfo.iscale = [xr, yr];
%       %iinfo.iscaleFormat = 'pq';
%       iinfo.iunit =  iinfo.imetadata.get('ResolutionUnit');
% 
%       return
%    end
% 
%    
%    %ZVI: todo: depth
%    xs = iinfo.imetadata.get('Scale Factor for Y');
%    if ~isempty(xs)
%       iinfo.iscale = [str2double(xs), str2double(iinfo.imetadata.get('Scale Factor for X'))];
%       iinfo.iunit  = 'um';
%    end
% 
%    %TODO:  %LIF
%    % ICS todo: check
%    voxelsizes = iinfo.imetadata.get('parameter scale');
%    if ~isempty(voxelsizes)% if possible pixelsizes are added (only ics files)
%       voxelsizes=str2num(voxelsizes); %#ok<ST2NM>
% 
%       if voxelsizes>1
%          iinfo.iscale = voxelsizes(2:5);
%          %iinfo.iscaleFormat = 'pqlt';
%       end
% 
%       return
%    end
  


% 
%    % TIF PhotometricInterpretation gives color interpretation
%    c = iinfo.imetadata.get('PhotometricInterpretation');
%    if ~isempty(c)
%       switch c
%          case 'BlackIsZero'
%             c = {'gray'};
%          case 'RGB'
%             c = {'r', 'g', 'b'};
%          otherwise
%             c = imcolorlist(iinfo.datasizeC);
%       end
% 
%       iinfo.icolor = c;
% 
%       return
%    end
% 
%    iinfo.icolor = imcolorlist(iinfo.datasizeC);
%    
%    %TODO:
%    %ZVI
%    %LIF


   
            
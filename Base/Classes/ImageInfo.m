classdef ImageInfo < matlab.mixin.Copyable
   
   % class that manages information about a single image in a file
   % for multiple series in a file choose ImageSourceList
   
   % notes:
   %    the info class is based on the info needed to read images via imreadBF / loci tools
   %    genreal organization:
   %                the ouput is represented in the data and cell size/format information
   %                raw input formats are stored in the raw data dn cell size / formats
   %                and possible reshaping of the raw cell data is done via the raw cell reshape options
   %                data flow: read image gives raw data when cell ranges are specified to be a singelton
   %                           the cell organizaton for reading is in the raw cell format
   %                           the cell data is then reshaped and reformatted into the output data and cells
   %                           with possible conversions between cell and data dimensions
   %                indexing / ranges:  ranges are index ranges or range keys
   
   properties
      % output data/cell sizes and formats when returned by the data and cell commands 
      idatasize      = [];           % size of the image 
      idataformat    = '';           % format of the image 

      idataclass = '';               % class of the image 
      
      icellsize   = 1;               % size of the cell structure
      icellformat = '';              % format of the cell structure
    
      
      % coordinate ranges
      irange = struct();             % specify a subset of coordinate ranges
      ikey   = struct();             % keys for certain coordinate ranges (e.g.  .C = {'GFP, 'DAPI'} corresponds to C = [1,2] etc)
      
 
      % input  data/cell sizes and formats when load from disk / stream etc
      irawdatasize = [];             % size of the raw data, [] = idatasize
      irawdataformat = '';           % format of raw data obtained by reading routines
      
      irawdataclass = '';            % class of the image 
      
      irawcellsize = [];             % size of the raw cell structure
      irawcellformat = '';           % format of the raw cells structure
      
      % reshaping
      ireshapefrom = {};             % formats dims in raw for reshaping 
      ireshapeto   = {};             % formats before and after reshaping each raw data array
      ireshapesize = {};             % sizes for the reshaping
           
      % meta
      ichannelname = [];             % channel name
      icolor = {};                   % colors for the color channel C/c
      icoloralpha = [];              % (optional) alpha values for the color channel C/c (=1 if no specified) 
  
      % intensity scale
      imaxintensity = [];            % maximal intensity of data
      iminintensity = [];            % minimal intensity of data
      
      imaxrawintensity = [];         % maximal intensity of raw data
      iminrawintensity = [];         % minimal intensity of raw data
      
      % spatial scales
      iscale = containers.Map;       % (optional) spatial scale of image in distance/time per pixel (coords: X, Y, Z, T)
      iunit  = containers.Map;       % (optional) spatial unit of image in distance/time per pixel (coords: X, Y, Z, T)
      
      iname  = [];                   % (optional) name of the image  
              
      imetadata = [];                % (optional) meta data
   end
   
   methods

      function obj = ImageInfo(varargin) % constructor
         %
         % ImageInfo()
         % ImageInfo(data)
         % ImageInfo(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageInfo') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.fromData(varargin{1});
            elseif ischar(varargin{1})
               obj.fromFile(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
         else
            obj = obj.fromParameter(varargin);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, 'i', varargin);
         obj.setCellDataSizeAndFormatToRaw();
      end
      
      function obj = fromData(obj, data)
         obj = imfrmtInfoFromData(data, obj);
         obj.setCellDataSizeAndFormatToRaw();
      end

      function obj = fromFile(obj, filename)
         info = imreadBFInfo(filename);
         obj.fromImageInfo(info);
         obj.setCellDataSizeAndFormatToRaw();
      end
      
      function obj = fromImageInfo(obj, info)
         for n = properties(ImageInfo)'
            obj.(n{1}) = info.(n{1});
         end
         obj.initializeRangeKeyFromChannelName;
      end
      
      function obj = fromImfInfo(obj, info)
         info = imfinfo2info(info);
         obj.fromImageInfo(info);
         obj.setCellDataSizeAndFormatToRaw();
      end

      function obj = setCellDataSizeAndFormatToRaw(obj)
         if isempty(obj.irawdataformat)
            obj.irawdataformat = obj.idataformat;
            obj.irawdatasize   = obj.idatasize;
         end
         
         if isempty(obj.irawcellformat)
            obj.irawcellformat = obj.icellformat;
            obj.irawcellsize   = obj.icellsize;
         end
      end
      
      function obj = setCellDataSizeAndFormatFromRaw(obj)
         if isempty(obj.idataformat)
            obj.idataformat = obj.irawdataformat;
            obj.idatasize   = obj.irawdatasize;
         end
         
         if isempty(obj.icellformat)
            obj.icellformat = obj.irawcellformat;
            obj.icellsize   = obj.irawcellsize;
         end
      end
      
      function obj = initializeDataClassFromRaw(obj)
         obj.idataclass = obj.irawdataclass;
      end
       
      function obj = initializeCellDataSizeFromRaw(obj)
         %
         % obj = initializeCellDataSizeFromRaw(obj)
         %
         % decription:
         %   sets data and cell size using raw data and cell size and format info 
         
 
         % reshaped sizes
         [obj.idatasize, obj.icellsize] = imfrmtReshapeCellDataSize( obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat, ...
                                                                     obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);    
                                                                  
         obj.irange = imfrmtReshapeRange(obj.rawCellDataSize, obj.rawCellDataFormat, obj.cellDataFormat, ...
                                         obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize, obj.irange);                                 
              
         % correct data sizes for ranges                                                             
         if ~isempty(obj.irange)
            obj.idatasize = imfrmtRangeSize(obj.dataSize, obj.dataFormat, obj.irange);
            obj.icellsize = imfrmtRangeSize(obj.cellSize, obj.cellFormat, obj.irange);
         end
      end
      
      function obj = initializeCellDataFormatFromRaw(obj)
         obj.idataformat = obj.irawdataformat;
         obj.icellformat = obj.irawcellformat;
      end
      
      function obj = initializeCellDataSizeAndFormatFromRaw(obj)
         obj.initializeCellDataFormatFromRaw();
         obj.initializeCellDataSizeFromRaw();
         obj.initializeDataClassFromRaw();
      end

      
      function obj = initializeReshape(obj)
         if ~iscell(obj.ireshapefrom)
            obj.ireshapefrom = {obj.ireshapefrom};
         end
         
         if ~iscell(obj.ireshapeto)
            obj.ireshapeto = {obj.ireshapeto};
         end
         
         if ~iscell(obj.ireshapesize)
            obj.ireshapesize = {obj.ireshapesize};
         end

         if length(obj.ireshapesize) ~= length(obj.ireshapefrom) || length(obj.ireshapesize) ~= length(obj.ireshapeto)
            error('ImageInfo: initializeReshape: inconsistent reshape formats or sizes');
         end
         
         [obj.idataformat, obj.icellformat] = imfrmtReshapeCellDataFormat(obj.irawdataformat, obj.irawcellformat, obj.ireshapefrom, obj.ireshapeto);

         obj.initializeCellDataSizeFromRaw;
         
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end
      
      
      function obj = initializeRangeKeyFromChannelName(obj)
         if ~isempty(obj.ichannelname)
            obj.ikey.C = obj.ichannelname;
         end
      end
      
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% keys
      
      
      function r = key(obj, varargin)
         if nargin == 1
            r = obj.ikey;
         else
            r = imfrmtRangeFromFormatAndVarargin(varargin{1}, obj.ikey);
         end
      end
      
      function obj = setKey(obj, varargin)
         obj.ikey = parseParameter(varargin);
      end
      
      function obj = addKey(obj, varargin)
         obj.ikey = parseParameter(obj.ikey, varargin);
      end
      
      function obj = resetKey(obj, varargin)
         obj.ikey = struct;
      end      
      

      function r = indexRangeFromKey(obj, varargin)
         % replaces keys with indices in any range (ouput and raw ranges)
         r = imfrmtRangeToIndexRange(obj.ikey, varargin);
      end
      
      function r = keyFromIndexRange(obj, varargin)
         % replaces indices with keys in any range (ouput and raw ranges)
         r = imfrmtRangeFromIndexRange(obj.ikey, varargin);
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% ranges 
      
      % note: irange and ikey are w.r.t. the full data format and size
      %       ranges restric data and cell size 

      
      function rn = rangeNames(obj)
         rn = fieldnames(obj.irange);
      end
      
      function obj = setRange(obj, varargin)
         obj.irange = struct;
         obj.irange = obj.rangeFromVarargin(varargin{:});
         obj.irange = imfrmtSortRange(obj.irange);
         
         obj.initializeCellDataSizeFromRaw;
         
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end
      
      function obj = addRange(obj, varargin)
         obj.irange = obj.rangeFromVarargin(varargin{:});
         obj.initializeCellDataSizeFromRaw;
         
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end
      
      function obj = resetRange(obj)
         obj.irange = struct;
         obj.initializeCellDataSizeFromRaw;
         
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end
      
      function obj = removeRange(obj, frmt)
         obj.irange = imfrmtRemoveRange(obj.irange, frmt);
         obj.initializeCellDataSizeFromRaw;
         
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end
      
      
      
      function range = range(obj, varargin)
         range = obj.rangeFromVarargin(varargin{:});
      end
      
      
      function [rRange, rReshape] = rawRange(obj, varargin)
         % determine raw range and optional reshaping sizes
         % inputs are converted via range keys
         % numerical indices are asusmed to be w.r.t. the output cell format / size
         % reshaping is w.r.t. to the full data format consistent with range specs
         
         range = obj.rangeFromVarargin(varargin{:});
          
         % inverse reshaping
         [rRange, rReshape] = imfrmtReshapeInverseCellDataRange(obj.fullDataSize, obj.fullDataFormat, obj.rawDataFormat, ...
                                                                obj.fullCellSize, obj.fullCellFormat, obj.rawCellFormat, ...
                                                                obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize, range);                                               
      end

      
      function range = rangeFromRawRange(obj, varargin)
         range = obj.rawRangeFromRawVarargin(varargin);
         range = imfrmtReshapeRange(obj.rawCellDataSize, obj.rawCellDataFormat, obj.cellDataFormat, ...
                                    obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize, range);     
      end

      
      function range = rawRangeFromRawRange(obj, varargin)
         range = obj.rawRangeFromRawVarargin(varargin);
      end
      
      
      % parsing routines
      
      function range = rangeFromCellIndex(obj, varargin)
         %
         % range = rangeFromCellIndex(obj)
         % range = rangeFromCellIndex(obj, idx)
         % range = rangeFromCellIndex(obj, sub1, sub2, ...)
         %
         % description:
         %     converts cell indices w.r.t. to the restricted cell size / format to a range
         %
         % input:
         %     idx    array of indices w.r.t to the restricted cell size and form
         %     sub*   sub array indices w.r.t. the restircted cell size and form 
         %
         % output:
         %     index range w.r.t. the full data cell  size / format
         
         % parse indices to restricted range
         range = imfrmtRangeFromIndex(obj.cellSize, obj.cellFormat, varargin{:});
         % convert restricted range to full index range
         range = imfrmtIndexRangeFromIndexRange(obj.irange, range);
         % complement to full range specs
         range = imfrmtRangeFromVarargin(obj.irange, range);
      end
      
      function range = rangeFromDataIndex(obj, varargin)
         %
         % range = rangeFromDataIndex(obj)
         % range = rangeFromDataIndex(obj, idx)
         % range = rangeFromDataIndex(obj, sub1, sub2, ...)
         %
         % description:
         %     converts data indices w.r.t. to the restricted data size / format to a range
         %
         % input:
         %     idx    array of indices w.r.t to the restricted data size and form
         %     sub*   sub array indices w.r.t. the restircted data size and form 
         %
         % output:
         %     index range w.r.t. the full data cell  size / format
         
         % parse indices to restricted range
         range = imfrmtRangeFromIndex(obj.dataSize, obj.dataFormat, varargin{:});
         % convert restricted range to full index range
         range = imfrmtIndexRangeFromIndexRange(obj.irange, range);
         % complement to full range specs
         range = imfrmtRangeFromVarargin(obj.irange, range);
      end
      
      
      function range = rangeFromVarargin(obj, varargin)
         %
         % range = rangeFromVarargin(obj, varargin)
         %
         % description:
         %     parses the range from inputs 
         %     numerical indices are w.r.t to the restricted cell format and size
         %     range specs are w.r.t. the full data / cell size / format
         %
         % output:
         %     index range w.r.t. the full data cell  size / format
         
         if nargin == 1 % nothing spcified return full specified range
            range = obj.irange;
 
         elseif nargin > 1 && isnumeric(varargin{1})  % form index 
            range = obj.rangeFromCellIndex(varargin{:});
            range = imfrmtSortRange(range);
         else
            range = imfrmtRangeFromVarargin(varargin); % from range specs
            range = imfrmtRangeToIndexRange(obj.key, range); 
            range = imfrmtRangeFromVarargin(obj.irange, range);
            range = imfrmtSortRange(range);
         end
      end
   
      
      function range = rawRangeFromRawCellIndex(obj, varargin)
         %
         % range = rawRangeFromRawCellIndex(obj)
         % range = rawRangeFromRawCellIndex(obj, idx)
         % range = rawRangeFromRawCellIndex(obj, sub1, sub2, ...)
         %
         % description:
         %     converts cell indices w.r.t. to the raw cell size / format to a raw range
         %
         % input:
         %     idx    array of indices w.r.t to the restricted raw cell size and form
         %     sub*   sub array indices w.r.t. the restircted raw cell size and form 
         %
         % output:
         %     range  raw index range w.r.t. the full data cell  size / format
         
         % parse indices to restricted range
         range = imfrmtRangeFromIndex(obj.rawCellSize, obj.rawCellFormat, varargin{:});

         % complement to full range specs
         range = imfrmtRangeFromVarargin(obj.rawRange, range);
      end
      
      
      function range = rawRangeFromRawVarargin(obj, varargin)
         % parses the raw range from inputs 
         % indices are w.r.t to raw cell format and size
         % all inputs are transformed via range keys
         
         if nargin == 1 % nothing spcified return full specified range
            range = obj.rawRange;

         elseif nargin > 1 && isnumeric(varargin{1})
            range = obj.rawRangeFromRawCellIndex(varargin{:});
            range = imfrmtSortRange(range);
         else
            range = imfrmtRangeFromVarargin(varargin); % from range specs 
            range = imfrmtRangeFromVarargin(obj.rawRange, range);
            range = imfrmtSortRange(range);
         end 
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data

      function d = dataDims(obj)
         d = length(obj.idataformat);
      end 

      function s = dataSize(obj, varargin) 
         if nargin == 1
            s = obj.idatasize;
         elseif nargin == 2 && ischar(varargin{1})
            s = imfrmtReformatSize(obj.idatasize, obj.idataformat, varargin{1});
         else
            s = imfrmtRangeSize(obj.dataSizeFromRawSize, obj.idataformat, obj.range(varargin{:}));
         end
      end

      function s = dataSizeFromRawSize(obj) 
         % returns full data size not constrained by ranges but using current data format
         s = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       obj.dataFormat,  obj.cellFormat,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function [s, frmt] = fullDataSize(obj) 
         % return full data size not constrained by ranges using full data format obtained from raw format
         frmt = obj.fullDataFormat;
         s = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       frmt,            obj.fullCellFormat,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
        
      function ds = dataSizeXYZCT(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'XYZCT');
      end

      function ds = dataSizeX(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'X');
      end
      function ds = dataSizeY(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'Y');
      end
      function ds = dataSizeZ(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'Z');
      end
      function ds = dataSizeC(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'C');
      end
      function ds = dataSizeT(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'T');
      end

      function df = dataFormat(obj)
         df = obj.idataformat;
      end
      
      function df = fullDataFormat(obj)
         % uses reshaping and raw format to get full version of the data frmt
         df = imfrmtReshapeFormat(obj.rawDataFormat, obj.reshapeFrom, obj.reshapeTo); 
         %df = imfrmtFlip(obj.dataFormat, df);
      end

      
      function dc = dataClass(obj)
         dc = obj.idataclass;
      end

%       function obj = setDataSize(obj, newsize)
%          obj.idatasize = newsize;
%       end
      
      function obj = setDataFormat(obj, newfrmt)
         %obj.idataformat = newfrmt;  
         obj = obj.reformatDataFormat(newfrmt);  % inconsistent bu more intuitive
      end
      
%       function obj = setDataFormatAndSize(obj, newfrmt, newsize)
%          obj.idatasize = newsize;
%          obj.idataformat = newfrmt;
%       end

      function obj = reformatDataFormat(obj, newfrmt)
         %obj.reformatColor(obj.idataformat, newfrmt);
         obj.idatasize   = imfrmtReformatSize(obj.fullDataSize, obj.fullDataFormat, newfrmt);
         obj.idataformat = newfrmt;
         
         obj.initializeCellDataSizeFromRaw;
         
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end

      function obj = renameDataFormat(obj, oldlab, newlab) 
          %
          % obj = renameDataFormat(obj, oldlab, newlab) 
          %
          % description:
          %    exchanges names of the dimensions oldlab(i) -> newlab(i)
          %    
          % input:
          %    oldlab     the format dimensions to change
          %    newlab     the new format dimensions 
           
          % check for conflict
          if any(ismember(setdiff(obj.idataformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameDataFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.idataformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.idataformat, oldlab);
          obj.idataformat(ids) = newlab;
          
          obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
          obj.ikey = renamestruct(obj.ikey, num2cell(oldlab), num2cell(newlab)); 
      end
       
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% cell

      function d = cellDim(obj)
         d = length(obj.icellformat);
      end

      function f = cellFormat(obj)
         f = obj.icellformat;
      end
      
      function s = cellSize(obj, varargin) 
         if nargin == 1
            s = obj.icellsize;
         elseif nargin == 2 && ischar(varargin{1})
            s = imfrmtReformatSize(obj.icellsize, obj.icellformat, varargin{1});
         else
            range = obj.rangeFromVarargin(varargin{:});
            s = imfrmtRangeSize(obj.cellSizeFromRawSize, obj.icellformat, range);
         end  
      end

      function s = cellSizeFromRawSize(obj) 
         % return full data size not constrained by ranges or curent data format
         [~, s] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       obj.dataFormat, obj.cellFormat,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end

      function [s, frmt] = fullCellSize(obj) 
         % return full data size not constrained by ranges but using current data format
         frmt = obj.fullCellFormat;
         [~, s] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       obj.fullDataFormat, frmt,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end

%       function obj = setCellSize(obj, newsize)
%          obj.isize = newsize;
%       end
         
      function df = fullCellFormat(obj)
         % uses reshaping and raw format info to get a version of the full data frmt
         df = imfrmtReshapeFormat(obj.rawCellFormat, obj.reshapeFrom, obj.reshapeTo); 
         %df = imfrmtFlip(obj.cellFormat, df);
      end
      
      function obj = setCellFormat(obj, newfrmt)
         %obj.icellformat = newfrmt;
         
         obj = obj.reformatCellFormat(newfrmt); %its inconsistent but more intuitive
      end

      function obj = reformatCellFormat(obj, newfrmt)
         obj.icellsize  = imfrmtReformatSize(obj.fullCellSize, obj.fullCellFormat, newfrmt);
         obj.icellformat = newfrmt;
         
         %obj.initializeCellDataSizeFromRaw;
         if isa(obj, 'ImageSource')
            obj.clearCellDataCache;
         end
      end 
      
      function obj = renameCellFormat(obj, oldlab, newlab)            
          % check for conflict
          if any(ismember(setdiff(obj.icellformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameCellFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.icellformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.icellformat, oldlab);
          obj.icellformat(ids) = newlab;
          
          obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
          obj.ikey = renamestruct(obj.ikey, num2cell(oldlab), num2cell(newlab)); 
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% raw data
        
 
      function d = rawDataDims(obj)
         d = length(obj.rawDataFormat);
      end 

      function rs = rawDataSize(obj, varargin)
         if nargin == 1
            rs = obj.irawdatasize;
         else
            rs = obj.rawRangeFromRawVarargin(varargin{:});
            rs = imfrmtRangeSize(obj.irawdatasize, obj.irawdataformat, rs);
         end
      end
      
      function rf = rawDataFormat(obj)
         rf = obj.irawdataformat;
      end
      
      function rc = rawDataClass(obj)
         rc = obj.irawdataclass;
      end
 
%       function obj = setRawDataSize(obj, newsize)
%          obj.irawdatasize = newsize;
%       end
      
      function obj = setRawDataFormat(obj, newfrmt)
         %obj.irawdataformat = newfrmt;
         obj.reformatRawDataFormat(newfrmt);
         
         
         if isa(obj, 'ImageSource')
            obj.clearRawCellDataCache;
         end
      end
      
      function obj = reformatRawDataFormat(obj, newfrmt)
         obj.irawdatasize   = imfrmtReformatSize(obj.rawDataSize, obj.rawDataFormat, newfrmt);
         obj.irawdataformat = newfrmt;
         
         
         if isa(obj, 'ImageSource')
            obj.clearRawCellDataCache;
         end
      end
      
      function obj = renameRawDataFormat(obj, oldlab, newlab)            
          % check for conflict
          if any(ismember(setdiff(obj.irawdataformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameRawDataFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.irawdataformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.irawdataformat, oldlab);
          obj.irawdataformat(ids) = newlab;
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% raw cell

      function d = rawCellDims(obj)
         d = length(obj.rawCellFormat);
      end 
      
      function rs = rawCellSize(obj, varargin)
         if nargin == 1
            rs = obj.irawcellsize;
         else
            rs = obj.rawRangeFromRawVarargin(varargin{:});
            rs = imfrmtRangeSize(obj.irawcellsize, obj.irawcellformat, rs);
         end
      end
      
      function rf = rawCellFormat(obj)
         rf = obj.irawcellformat;
      end
      
%       function obj = setRawCellSize(obj, newsize)
%          obj.irawcellsize = newsize;
%       end
      
      function obj = setRawCellFormat(obj, newfrmt)
         %obj.irawcellformat = newfrmt;
         obj.reformatRawCellFormat(newfrmt); % consistent with intuitive setCellFormat
       
         
         if isa(obj, 'ImageSource')
            obj.clearRawCellDataCache;
         end
      end
      
      
      function obj = reformatRawCellFormat(obj, newfrmt)
         obj.irawcellsize  = imfrmtReformatSize(obj.rawCellSize, obj.rawCellFormat, newfrmt);
         obj.irawcellformat = newfrmt;
         
         
         if isa(obj, 'ImageSource')
            obj.clearRawCellDataCache;
         end
      end
      
      function obj = renameRawCellFormat(obj, oldlab, newlab)            
          % check for conflict
          if any(ismember(setdiff(obj.irawcellformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameRawCellFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.irawcellformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.irawcellformat, oldlab);
          obj.irawcellformat(ids) = newlab;
      end
      

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% cell data 
      
      function f = cellDataFormat(obj)
         f = [obj.dataFormat, obj.cellFormat];
      end

      function f = cellDataSize(obj, varargin)
         f = [obj.dataSize(varargin{:}), obj.cellSize(varargin{:})];
      end
          
      function f =rawCellDataFormat(obj)
         f = [obj.rawDataFormat, obj.rawCellFormat];
      end

      function f = rawCellDataSize(obj, varargin)
         f = [obj.rawDataSize(varargin{:}), obj.rawCellSize(varargin{:})];
      end


      function [ds, cs] = cellDataSizeFromRaw(obj) 
         % return full cell / data size not constrained by ranges but using current data format
         [ds, cs] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                             obj.dataFormat, obj.cellFormat,...
                                             obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end

      function [ds, cs, df, cf] = fullCellDataSize(obj) 
         % return full cell / data size not constrained by ranges or curent data format
         [df, cf] = imfrmtReshapeCellDataFormat(obj.rawDataFormat, obj.rawCellFormat, obj.reshapeFrom. obj.reshapeTo);  
         [ds, cs] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       df, cf,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function [df, cf] = fullCellDataFormat(obj)
         % uses reshaping and raw format info to get a version of the full data frmt
         [df, cf] = imfrmtReshapeCellDataFormat(obj.rawDataFormat, obj.rawCellFormat, obj.reshapeFrom. obj.reshapeTo); 
         df = imfrmtFlip(obj.dataFormat, df);
         cf = imfrmtFlip(obj.cellFormat, cf);
      end
      

      function obj = setCellDataFormat(obj, newdatafrmt, newcellfrmt)    
         obj.setDataFormat(newdatafrmt);
         obj.setCellFormat(newcellfrmt);
      end
 
      function obj = renameCellDataFormat(obj, oldlab, newlab)
         obj.renameDataFormat(oldlab, newlab);
         obj.renameCellFormat(oldlab, newlab);
      end
      
      function obj = renameRawCellDataFormat(obj, oldlab, newlab)
         obj.renameRawDataFormat(oldlab, newlab);
         obj.renameRawCellFormat(oldlab, newlab);
      end
      
      
%       function obj = renameFormat(obj, oldlab, newlab)
%          obj.renameDataFormat(oldlab, newlab);
%          obj.renameCellFormat(oldlab, newlab);
%          
%          obj.renameRawDataFormat(oldlab, newlab);
%          obj.renameRawCellFormat(oldlab, newlab);
%          
%          % note: no renamming of reshape formats done here
%       end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% indices
      
      % note: three sorts of indices: 
      %          - out format indices w.r.t to the restricted range
      %          - full out format indices not restricted by range
      %          - raw indices
      %
      %          - keys and internal ranges are alsways w.r.t the full range

      
      % data 
            
      function n = nDataIndices(obj, varargin)
         n = prod(obj.dataSize(obj, varargin));
      end
      
      function n = nFullDataIndices(obj, varargin)
         n = prod(obj.fullDataSize(varargin{:}));
      end
      
            
      function n = nCellIndices(obj, varargin)
         n = prod(obj.cellSize(varargin{:}));
      end

      function n = nFullCellIndices(obj, varargin)
         n = prod(obj.fullCellSize(varargin{:}));
      end
      
      function n = nCells(obj, varargin)
         n = prod(obj.cellSize(varargin{:}));
      end
      
      
      function id = dataIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            id = imfrmtIndexFromSizeAndVarargin(obj.dataSize, varargin{:});
         else
            range = obj.rangeFromVarargin(varargin{:});
            id = imfrmtRangeToIndex(obj.fullDataSize, obj.fullDataFormat, range);
            id = imfrmtFullIndexToIndex(obj.fullDataSize, obj.fullDataFormat, obj.irange, id); 
         end
      end

      function id = fullDataIndex(obj, varargin)
         % converts index w.r.t. to restricted range to index w.r.t. to full range 
         if nargin >= 2 && isnumeric(varargin{1})
            id = imfrmtIndexFromSizeAndVarargin(obj.fullDataSize, varargin{:});
         else
            range = obj.rangeFromVarargin(varargin{:});
            id = imfrmtRangeToIndex(obj.fullDataSize, obj.fullDataFormat, range);
         end
      end
      

      
      function id = cellIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            id = imfrmtIndexFromSizeAndVarargin(obj.cellSize, varargin{:});
         else
            range = obj.rangeFromVarargin(varargin{:});
            id = imfrmtRangeToIndex(obj.fullCellSize, obj.fullCellFormat, range);
            id = imfrmtFullIndexToIndex(obj.fullCellSize, obj.fullCellFormat, obj.irange, id);
         end
      end

      function id = fullCellIndex(obj, varargin)
         % converts index w.r.t. to restricted range to index w.r.t. to full range 
         if nargin >= 2 && isnumeric(varargin{1})
            id = imfrmtIndexFromSizeAndVarargin(obj.fullCellSize, varargin{:});
         else
            range = obj.rangeFromVarargin(varargin{:});
            id = imfrmtRangeToIndex(obj.fullCellSize, obj.fullCellFormat, range);
         end
      end
      
      
      
      function id = fullDataIndexToDataIndex(obj, id)
         id = imfrmtFullIndexToIndex(obj.fullDataSize, obj.fullDataFormat, obj.irange, id);
      end
      
      function id = fullCellIndexToCllIndex(obj, id)
         id = imfrmtFullIndexToIndex(obj.fullCellSize, obj.fullCellFormat, obj.irange, id);
      end
      
      
      
      function id = dataIndexToSubIndex(obj, id, varargin)
         if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
            id = [id, varargin{:}];
         else
            id = imind2sub(obj.dataSize, id);
         end
      end

      function id = cellIndexToSubIndex(obj, id, varargin)
         if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
            id = [id, varargin{:}];
         else
            id = imind2sub(obj.cellSize, id);
         end
      end
      
                
%       function range = dataIndexToRange(obj, varargin)
%          % returns ranges restricted by interal range settings and data specs in varargin
%          range = imfrmtRangeFromIndex(obj.dataSize, obj.dataFormat, varargin{:});
%          range = imfrmtRangeFromVarargin(obj.irange, range);
%       end



      function n = nRawDataIndices(obj, varargin)
         n = prod(obj.rawDataSize(varargin{:}));
      end
      
      function n = nRawCellIndices(obj, varargin)
         n = prod(obj.rawCellSize(varargin{:}));
      end

      function n = nRawCells(obj, varargin)
         n = prod(obj.rawCellSize(varargin{:}));
      end
      
  
      function id = rawDataIndexFromRawVarargin(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            id = imfrmtIndexFromVarargin(obj.rawDataSize, varargin{:});
         else
            range = obj.rawRangeFromRawVarargin(varargin{:});
            id = imfrmtRangeToIndex(obj.rawDataSize, obj.rawDataFormat, range);
         end
      end

      function id = rawCellIndexFromRawVarargin(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            id = imfrmtIndexFromVarargin(obj.rawCellSize, varargin{:});
         else
            range = obj.rawRangeFromRawVarargin(varargin{:});
            id = imfrmtRangeToIndex(obj.rawCellSize, obj.rawCellFormat, range);
         end
      end
      
      
      function id = rawDataIndex(obj, varargin)
         range = obj.rawRange(varargin{:})
         id = imfrmtRangeToIndex(obj.rawDataSize, obj.rawDataFormat, range);
      end

      function id = rawCellIndex(obj, varargin)
         range = obj.rawRange(varargin{:});
         id = imfrmtRangeToIndex(obj.rawCellSize, obj.rawCellFormat, range);
      end
      
      function range = rawDataIndexToRange(obj, varargin)
         range = imfrmtRangeFromIndex(obj.rawDataSize, obj.rawDataFormat, varargin{:});
         range = imfrmtRangeFromVarargin(obj.rawRange, range);
      end
      
      function range = rawCellIndexToRange(obj, varargin)
         range = imfrmtRangeFromIndex(obj.rawCellSize, obj.rawCellFormat, varargin{:});
         range = imfrmtRangeFromVarargin(obj.rawRange, range);
      end
     
      
  
%       function coords = dataIndexToCoordinate(obj, id, varargin)
%          if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
%             id = imsub2ind(obj.dataSize, [id, varargin{:}]);
%          end 
%          coords = imfrmtIndexToCoordinate(obj.dataSize, obj.dataFormat, id);
%       end
%       
%       function coords = cellIndexToCoordinate(obj, id, varargin)
%          if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
%             id = imsub2ind(obj.cellSize, [id, varargin{:}]);
%          end 
%          coords = imfrmtIndexToCoordinate(obj.cellSize, obj.cellFormat, id);
%       end
%       
%       
%       function coords = rawDataIndexToCoordinate(obj, id, varargin)
%          if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
%             id = imsub2ind(obj.rawDataSize, [id, varargin{:}]);
%          end 
%          coords = imfrmtIndexToCoordinate(obj.rawDataSize, obj.rawDataFormat, id);
%       end
%       
%       function coords = rawCellIndexToCoordinate(obj, id, varargin)
%          if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
%             id = imsub2ind(obj.rawCellSize, [id, varargin{:}]);
%          end 
%          coords = imfrmtIndexToCoordinate(obj.rawCellSize, obj.rawCellFormat, id);
%       end
       
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% reshaping 
          
      function r = reshapeFrom(obj)
         r = obj.ireshapefrom;
      end
      
      function r = reshapeTo(obj)
         r = obj.ireshapeto;
      end
      
      function r = reshapeSize(obj)
         r = obj.ireshapesize;
      end     
      
      function obj = setReshape(obj, reshapeFrom, reshapeTo, reshapeSize)
         obj.ireshapefrom = reshapeFrom;
         obj.ireshapeto   = reshapeTo;
         obj.ireshapesize = reshapeSize;   
         obj.initializeReshape();
      end

      function obj = addReshape(obj, reshapeFrom, reshapeTo, reshapeSize)
         obj.ireshapefrom{end+1} = reshapeFrom;
         obj.ireshapeto{end+1}   = reshapeTo;
         obj.ireshapeSize{end+1} = reshapeSize;
         obj.initializeReshape();
      end
 
      function d = reshapeData(obj, d)
         % uses reshaping and format info to reshape the raw data d to the return data
         d = imfrmtReshape(d, obj.rawDataFormat, obj.dataFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);  
      end
      
      function d = reshapeCell(obj, d)
         % uses reshaping and format info to reshape the raw data d to the return data
         d = imfrmtReshape(d, obj.rawCellFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);  
      end

      function d = reshapeCellData(obj, d)
         % uses reshaping and format info to reshape the raw data d to the return data
         d = imfrmtReshapeCellData(d,obj.rawDataFormat,  obj.rawCellFormat, obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);  
      end
    
      


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% assignments
      
      function asgn = dataAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.dataSize, obj.dataFormat, varargin); 
      end
      
      function asgn = cellAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.cellSize, obj.cellFormat, varargin);
      end
      
      function asgn = rawDataAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.rawDataSize, obj.rawDataFormat, varargin); 
      end
      
      function asgn = rawCellAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.rawCellSize, obj.rawCellFormat, varargin);
      end
      
      
      function i = dataFormatPosition(obj, frmt)
         %
         % i = dataFormatPosition(obj, frmt)
         % 
         % description: 
         %    returns position of the format names in frmt in obj.dataFormat
         
         [~, i] = ismember(frmt, obj.dataFormat);
      end
      
      function i = cellFormatPosition(obj, frmt)
         [~, i] = ismember(frmt, obj.cellFormat);
      end
 
      function i = rawDataFormatPosition(obj, frmt)
         [~, i] = ismember(frmt, obj.rawDataFormat);
      end
      
      function i = rawCellFormatPosition(obj, frmt)
         [~, i] = ismember(frmt, obj.rawCellFormat);
      end
  
   
  

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% intensities

      function m = maxDataIntensity(obj, varargin)
         m = obj.imaxintensity;
      end
      
      function m = minDataIntensity(obj, varargin)
         m = obj.iminintensity;
      end
      
      function m = maxRawIntensity(obj, varargin)
         m = obj.imaxrawintensity;
      end
      
      function m = minRawIntensity(obj, varargin)
         m = obj.iminrawintensity;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% casting
      
      function d = castToDataClass(obj, d)
         if iscell(d)
            d = cellfunc(@(x) cast(x,obj.dataclass), d);
         else
            d = cast(d, obj.idataclass);
         end
      end
          
      function d = castToRawClass(obj, d)
         if iscell(d)
            d = cellfunc(@(x) cast(x,obj.rawdataclass), d);
         else
            d = cast(d, obj.irawdataclass);
         end
      end
      
      
      function d = castRawToData(obj, d)
         param.data.min = obj.minRawIntensity;
         param.data.max = obj.maxRawIntensity;
         param.rescale.max = obj.maxDataIntensity;
         param.rescale.min = obj.minDataIntenisty;
         param.class = obj.dataClass;
         
         if iscell(d)
            d = cellfunc(@(x) imrescale(x, param), d);
         else
            d = imrescale(d, param);
         end
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% scaling
      
      function sc = scale(obj)
         dfrmt = 'XYZT';
         pos = imfrmtPosition(dfrmt, obj.dataFormat);
         
         sc = zeros(1,length(pos));
         for i = 1
            sc(i) = obj.iscale(dfrmt(pos(i)));
         end
      end
      
      function sc = scaleX(obj)
         sc = obj.iscale('X');
      end
      
      function sc = scaleY(obj)
         sc = obj.iscale('Y');
      end
 
      function sc = scaleZ(obj)
         sc = obj.iscale('Z');
      end

      function sc = scaleT(obj)
         sc = obj.iscale('T');
      end
      
      function obj = setScale(obj, dim, scale)
         obj.iscale(upper(dim)) = scale;
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% units
      
      function sc = units(obj)
         dfrmt = 'XYZT';
         pos = imfrmtPosition(dfrmt, obj.dataFormat);
         
         sc = zeros(1,length(pos));
         for i = 1
            sc(i) = obj.iunit(dfrmt(pos(i)));
         end
      end
      
      function sc = unitX(obj)
         sc = obj.iunit('X');
      end
      
      function sc = unitY(obj)
         sc = obj.iunit('Y');
      end
 
      function sc = unitZ(obj)
         sc = obj.iunit('Z');
      end

      function sc = unitT(obj)
         sc = obj.iunit('T');
      end
     
      function obj = setUnit(obj, dim, unit)
         obj.iunit(upper(dim)) = unit;
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% names
      
      function n = name(obj)
         n = obj.iname;
      end
      
      function obj = setName(obj, name)
         obj.iname = name;
      end
          
      function n = channelName(obj, varargin)
         n = obj.ichannelname(varargin{:});
      end
      
      function obj = setChannelName(obj, varargin)
         if nargin == 3
            obj.ichannelname(varargin{1}) = varargin{2};
         else
            obj.ichannelname = varargin{1};
         end
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% meta
      
      function md = metaData(obj, varargin)
         if nargin > 1
            md = obj.imetadata;
            i = find(ismember(md.Parameter, varargin{1}));
            if ~isempty(i)
               md = md.Value(i);
            end
         else
            md = obj.imetadata;
         end
      end 
      
      
      function mdn = metaDataNames(obj)
         mdn = obj.imetadata;
         mdn = mdn.Parameter;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% colors
      
      function n = color(obj, varargin)
         n = obj.icolor(varargin{:});
      end
      
      function obj = setColor(obj, varargin)
          if nargin == 3
            obj.icolor(varargin{1}) = varargin{2};
          else
            if ischar(varargin{1})
               obj.icolor = varargin(1);
            else
               obj.icolor = varargin{1};
            end
         end
      end
      
            
      function n = colorAlpha(obj, varargin)
         n = obj.icoloralpha(varargin{:});
      end
      
      function obj = setColorAlpha(obj, varargin)
         if nargin == 3
            obj.icoloralpha(varargin{1}) = varargin{2};
         else
            obj.icoloralpha = varargin{1};
         end
      end

      function obj = reformatColors(obj, ifrmt, ofrmt)
         i = find(lower(ifrmt)== 'c', 1, 'first');
         if ~isempty(i)
            o = find(lower(ofrmt)== 'c', 1, 'first');
            if ~isempty(o)
               if ifrmt(i) ~= ofrmt(o)
                  if ~isempty(obj.icolor)
                     obj.icolor = flip(obj.icolor);
                  end
                  if ~isempty(obj.icoloralpha)
                     obj.ialphacolor = flip(obj.icoloralpha);
                  end
                  if ~isempty(obj.ichannelnames)
                     obj.ichannelnames = flip(obj.ichannelnames);
                  end
               end
            end
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% info

      function istr = infoString(obj, varargin)
         if nargin > 1
            cls = varargin{1};
         else
            cls = '';
         end
         istr = ['ImageSource: ',  cls];
         if ~isempty(obj.name)
            istr = [istr, '\nname:   ', var2char(obj.name)];
         end
         istr = [istr, '\nsize:       ', var2char(obj.dataSize)];
         istr = [istr, '\nformat:     ', var2char(obj.dataFormat)];
         istr = [istr, '\ncellsize    ', var2char(obj.cellSize)]; 
         istr = [istr, '\ncellformat: ', var2char(obj.cellFormat)];     
         istr = [istr, '\nclass:      ', var2char(obj.dataClass)];
         istr = [istr, '\ncolor:      ', var2char(obj.color)];
      end

      function printInfo(obj)
         fprintf([obj.infoString, '\n']);
      end
      
      
   end
end
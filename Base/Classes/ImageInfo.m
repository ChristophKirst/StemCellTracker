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
   
   properties
      % output data/cell sizes and formats when returned by the data and cell commands 
      idatasize      = [];           % size of the image 
      idataformat    = '';           % format of the image 

      idataclass = '';               % class of the image 
      
      icellsize   = 1;               % size of the cell structure
      icellformat = '';              % format of the cell structure
    
      
      % coordinate ranges
      irange = struct();             % specify a subset of coordinate ranges
      irangekey = struct();          % keys for certain coordinate ranges (e.g.  .C = {'GFP, 'DAPI'} corresponds to C = [1,2] etc)
      
 
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
      imaxintensity = [];            % minimal intensity of data
      iminintensity = [];            % maximal intensity of data
      
      imaxrawintensity = [];         % minimal intensity of raw data
      iminrawintensity = [];         % maximal intensity of raw data
      
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
         obj.initializeRawFormatsAndSizes();
      end
      
      function obj = fromData(obj, data)
         obj = imfrmtInfoFromData(data, obj);
         obj.initializeRawCellDataFormatsAndSizesFromData();
      end

      function obj = fromFile(obj, filename)
         info = imreadBFInfo(filename);
         obj.fromImageInfo(info);
         obj.initializeRawCellDataFormatsAndSizesFromData();
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
         obj.initializeRawCellDataFormatsAndSizesFromData();
      end

      function obj = initializeRawCellDataFormatsAndSizesFromData(obj)
         if isempty(obj.irawdataformat)
            obj.irawdataformat = obj.idataformat;
            obj.irawdatasize   = obj.idatasize;
         end
         
         if isempty(obj.irawcellformat)
            obj.irawcellformat = obj.icellformat;
            obj.irawcellsize   = obj.icellsize;
         end
      end
      
      function obj = initializeCellDataFormatsAndSizesFromRaw(obj)
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
       
      function obj = initializeDataAndCellSizeFromRaw(obj)
         %
         % obj = initializeDataAndCellSizeFromRaw(obj)
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
      
      function obj = initializeDataAndCellFormatFromRaw(obj)
         obj.idataformat = obj.irawdataformat;
         obj.icellformat = obj.irawcellformat;
      end
      
      function obj = initializeDataAndCellFromRaw(obj)
         obj.initializeDataAndCellFormatFromRaw();
         obj.initializeDataAndCellSizeFromRaw();
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

         obj.initializeDataAndCellSizeFromRaw;
         obj.clearCellDataCache;
      end
      
      
      function obj = initializeRangeKeyFromChannelName(obj)
         if ~isempty(obj.ichannelname)
            obj.irangekey.C = obj.ichannelname;
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
         %elseif nargin == 2 && ischar(varargin{1})
         %   s = imfrmtReformatSize(obj.idatasize, obj.idataformat, varargin{1});
         else
            range = obj.dataRangeFromVarargin(varargin{:});
            s = imfrmtRangeSize(obj.dataSizeFromRaw, obj.idataformat, range);
         end
      end

      function s = dataSizeFromRaw(obj) 
         % return full data size not constrained by ranges ubut using current data format
         s = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       obj.dataFormat,  obj.cellFormat,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function [s, frmt] = dataSizeFromRawFormat(obj) 
         % return full data size not constrained by ranges using full data format obtained from raw format
         frmt = obj.dataFormatFromRawFormat;
         s = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       frmt,            obj.cellFormatFromRawFormat,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function df = dataFormatFromRawFormat(obj)
         % uses reshaping and raw format to get full version of the data frmt
         df = imfrmtReshapeFormat(obj.rawDataFormat, obj.reshapeFrom. obj.reshapeTo); 
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
      
      function dc = dataClass(obj)
         dc = obj.idataclass;
      end
      
      
      function obj = setDataSize(obj, newsize)
         obj.idatasize = newsize;
      end
      
      function obj = setDataFormat(obj, newfrmt)
         obj.idataformat = newfrmt;  
      end
      
      function obj = setDataFormatAndSize(obj, newfrmt, newsize)
         obj.idatasize = newsize;
         obj.idataformat = newfrmt;
      end

      function obj = reformatDataFormat(obj, newfrmt)
         %obj.reformatColor(obj.idataformat, newfrmt);
         obj.idatasize   = imfrmtReformatSize(obj.dataSizeFromRaw, obj.idataformat, newfrmt);
         obj.idataformat = newfrmt;
         
         obj.initializeDataAndCellSizeFromRaw;
         obj.clearCellDataCache;
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
          obj.irangename = renamestruct(obj.irangename, num2cell(oldlab), num2cell(newlab)); 
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
         %elseif nargin == 2 && ischar(varargin{1})
         %   s = imfrmtReformatSize(obj.icellsize, obj.icellformat, varargin{1});
         else
            range = obj.cellRangeFromVarargin(varargin{:});
            s = imfrmtRangeSize(obj.cellSizeFromRaw, obj.icellformat, range);
         end  
      end

      function s = cellSizeFromRaw(obj) 
         % return full data size not constrained by ranges or curent data format
         [~, s] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       obj.dataFormat, obj.cellFormat,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end

      function [s, frmt] = cellSizeFromRawFormat(obj) 
         % return full data size not constrained by ranges but using current data format
         frmt = obj.cellFormatFromRawFormat;
         [~, s] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       obj.dataFormatFromRawFormat, frmt,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function df = cellFormatFromRawFormat(obj)
         % uses reshaping and raw format info to get a version of the full data frmt
         df = imfrmtReshapeFormat(obj.rawCellFormat, obj.reshapeFrom. obj.reshapeTo); 
      end
      
      
      function obj = setCellSize(obj, newsize)
         obj.isize = newsize;
      end
         
      function obj = setCellFormat(obj, newfrmt)
         obj.icellformat = newfrmt;
      end

      function obj = reformatCellFormat(obj, newfrmt)
         obj.icellsize  = imfrmtReformatSize(obj.cellSize, obj.cellFormat, newfrmt);
         obj.icellformat = newfrmt;
         
         obj.initializeDataAndCellSizeFromRaw;
         obj.clearCellDataCache;
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
            rs = obj.rawDataRangeFromVarargin(varargin{:});
            rs = imfrmtRangeSize(obj.irawdatasize, obj.irawdataformat, rs);
         end
      end
      
      function rf = rawDataFormat(obj)
         rf = obj.irawdataformat;
      end
      
      function rc = rawDataClass(obj)
         rc = obj.irawdataclass;
      end
 
      function obj = setRawDataSize(obj, newsize)
         obj.irawdatasize = newsize;
      end
      
      function obj = setRawDataFormat(obj, newfrmt)
         obj.irawdataformat = newfrmt;
      end
      
      function obj = reformatRawDataFormat(obj, newfrmt)
         obj.irawdatasize   = imfrmtReformatSize(obj.rawDataSize, obj.rawDataFormat, newfrmt);
         obj.irawdataformat = newfrmt;
         
         obj.clearRawCellDataCache;
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
            rs = obj.rawCellRangeFromVarargin(varargin{:});
            rs = imfrmtRangeSize(obj.irawcellsize, obj.irawcellformat, rs);
         end
      end
      
      function rf = rawCellFormat(obj)
         rf = obj.irawcellformat;
      end
      
      function obj = setRawCellSize(obj, newsize)
         obj.irawcellsize = newsize;
      end
      
      function obj = setRawCellFormat(obj, newfrmt)
         obj.irawcellformat = newfrmt;
      end
      
      function obj = reformatRawCellFormat(obj, newfrmt)
         obj.irawcellsize  = imfrmtReformatSize(obj.rawCellSize, obj.rawCellFormat, newfrmt);
         obj.irawcellformat = newfrmt;
         
         obj.clearRawCellDataCache;
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
      %%% cell data / full format 
      
      function f = cellDataFormat(obj)
         f = [obj.dataFormat, obj.cellFormat];
      end

      function f = cellDataSize(obj, varargin)
         f = [obj.dataSize(varargin{:}), obj.dataFormat(varargin{:})];
      end
          
      function f =rawCellDataFormat(obj)
         f = [obj.rawDataFormat, obj.rawCellFormat];
      end

      function f = rawCellDataSize(obj, varargin)
         f = [obj.rawDataSize(varargin{:}), obj.rawDataFormat(varargin{:})];
      end


      function [s, cs] = cellDataSizeFromRaw(obj) 
         % return full cell / data size not constrained by ranges but using current data format
         [s, cs] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                             obj.dataFormat, obj.cellFormat,...
                                             obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end

      function [ds, cs, df, cf] = cellDataSizeFromRawFormat(obj) 
         % return full cell / data size not constrained by ranges or curent data format
         [df, cf] = imfrmtReshapeCellDataFormat(obj.rawDataFormat, obj.rawCellFormat, obj.reshapeFrom. obj.reshapeTo);  
         [ds, cs] = imfrmtReshapeCellDataSize(obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat,...
                                       df, cf,...
                                       obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function [df, cf] = cellDataFormatFromRawFormat(obj)
         % uses reshaping and raw format info to get a version of the full data frmt
         [df, cf] = imfrmtReshapeCellDataFormat(obj.rawDataFormat, obj.rawCellFormat, obj.reshapeFrom. obj.reshapeTo); 
      end
      

      function obj = setCellDataFormat(obj, newdatafrmt, newcellfrmt)
         obj.idataformat = newdatafrmt;
         obj.icellformat = newcellfrmt;
         obj.initializeDataAndCellSizeFromRaw();
      end
      
       
      function obj = renameCellDataFormat(obj, oldlab, newlab)
         obj.renameDataFormat(oldlab, newlab);
         %obj.renameRawDataFormat(oldlab, newlab);
        
         obj.renameCellFormat(oldlab, newlab);
         %obj.renameRawCellFormat(oldlab, newlab);
         
         %obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
      end
      
      function obj = renameRawCellDataFormat(obj, oldlab, newlab)
         %obj.renameDataFormat(oldlab, newlab);
         obj.renameRawDataFormat(oldlab, newlab);
        
         %obj.renameCellFormat(oldlab, newlab);
         obj.renameRawCellFormat(oldlab, newlab);
         
         %obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
      end
      
      function obj = renameFormat(obj, oldlab, newlab)
         obj.renameDataFormat(oldlab, newlab);
         obj.renameRawDataFormat(oldlab, newlab);
        
         obj.renameCellFormat(oldlab, newlab);
         obj.renameRawCellFormat(oldlab, newlab);
         
         %obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
      end

      
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
      %%% ranges
          
      function r = range(obj, varargin)
         
         range = obj.cellRangeFromVarargin(varargin{:});    
         r = imfrmtParseRangeLast(obj.irange, range);
      end
 
      function [rRange, rReshape] = rawRange(obj, varargin)
         % determine raw range and optional reshaping sizes
         % inputs are converted via range keys
         % numerical indices are asusmed to be w.r.t. the cell format  / size
         
         range = obj.cellRangeFromVarargin(varargin{:});
         
         % inverse reshaping
         [rRange, rReshape] = imfrmtReshapeInverseCellDataRange(obj.dataSizeFromRaw, obj.dataFormat, obj.rawDataFormat, ...
                                                                obj.cellSizeFromRaw, obj.cellFormat, obj.rawCellFormat, obj.reshapeFrom, ...
                                                                obj.reshapeTo, obj.reshapeSize, range);                                               
      end
      
      function r = rangeKey(obj)
         r = obj.irangekey;
      end
      
      function r = rangeToIndexRange(obj, varargin)
         r = imfrmtRangeToIndexRange(obj.irangekey, parseParameter(varargin));
      end
      
      function r = rangeFromIndexRange(obj, varargin)
         r = imfrmtRangeToIndexRange(obj.irangekey, parseParameter(varargin));
      end
      
      function obj = setRangeKey(obj, varargin)
         obj.irangekey = parseParameter(varargin);
      end
      
      function obj = addRangeKey(obj, varargin)
         obj.irangekey = parseParameter(obj.irangekey, varargin);
      end
      
      function obj = resetRangeKey(obj, varargin)
         obj.irangekey = struct;
      end
      
 
      function rn = rangeNames(obj)
         rn = fieldnames(obj.irange);
      end
      
      function obj = setRange(obj, varargin)
         obj.irange = imfrmtParseRange(obj.cellDataSize, obj.cellDataFormat, varargin);
         obj.initializeDataAndCellSizeFromRaw;
         obj.clearCellDataCache;
      end
      
      function obj = addRange(obj, varargin)
         obj.irange = parseParameter(obj.irange, varargin);
         obj.initializeDataAndCellSizeFromRaw;
         obj.clearCellDataCache;
      end
      
      function obj = resetRange(obj)
         obj.irange = struct();
         obj.initializeDataAndCellSizeFromRaw;
         obj.clearCellDataCache;
      end
          
       
      function range = dataRangeFromVarargin(obj, varargin)
         % parses the raw range from inputs 
         % indices are w.r.t to the data format and size
         % all inputs are transformed via range keys
         
         if nargin > 1 && isnumeric(varargin{1})
            range = obj.dataIndexToCoordinate(varargin{:});
         else
            range = imfrmtParseRangeLast(obj.irange, obj.rangeToIndexRange(varargin));
         end 
      end
      
      function range = cellRangeFromVarargin(obj, varargin)
         % parses the raw range from inputs 
         % indices are w.r.t to the cell format and size
         % all inputs are transformed via range keys
         
         if nargin > 1 && isnumeric(varargin{1})
            range = obj.cellIndexToCoordinate(varargin{:});
            range = imfrmtParseRangeLast(obj.irange, range);
         else
            range = imfrmtParseRangeLast(obj.irange, obj.rangeToIndexRange(varargin));
         end 
      end
      
      function range = rawDataRangeFromVarargin(obj, varargin)
         % parses the raw range from inputs 
         % indices are w.r.t to the data format and size
         % all inputs are assumed to be index ranges (no transform via keys!)
         
         if nargin > 1 && isnumeric(varargin{1})
            range = obj.rawDataIndexToCoordinate(varargin{:});
         else
            %range = obj.rangeToIndexRange(varargin);
            %range = imfrmtParseRangeLast(obj.rawRange, range);
            range = imfrmtParseRangeLast(obj.rawRange, varargin{:});
         end 
      end
      
      function range = rawCellRangeFromVarargin(obj, varargin)
         % parses the raw range from inputs 
         % indices are w.r.t to the cell format and size
         % all inputs are assumed to be index ranges (no transform via keys!)

         if nargin > 1 && isnumeric(varargin{1})
            range = obj.rawCellIndexToCoordinate(varargin{:});
         else
            %range = obj.rangeToIndexRange(varargin);
            %range = imfrmtParseRangeLast(obj.rawRange, range);
            range = imfrmtParseRangeLast(obj.rawRange, varargin{:});
         end 
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
      %%% indices
      
            
      function n = nDataIndices(obj, varargin)
         n = prod(obj.dataSize(varargin{:}));
      end
      
      function n = nCellIndices(obj, varargin)
         n = prod(obj.cellSize(varargin{:}));
      end
      
      function n = nCells(obj, varargin)
         n = prod(obj.cellSize(varargin{:}));
      end
      
      function n = nRawDataIndices(obj, varargin)
         n = prod(obj.rawDataSize(varargin{:}));
      end
      
      function n = nRawCellIndices(obj, varargin)
         n = prod(obj.rawCellSize(varargin{:}));
      end

      function i = dataIndex(obj, varargin)
         % returns indices of the data corresponding to the input specs
         % ranges are converted via keys
         % index is wrt to the data size
         
         if nargin >= 2 && isnumeric(varargin{1}) && isscalar(varargin{1})
            if nargin > 2
               i = imsub2ind(obj.dataSize, [varargin{:}]);
            else
               i = varargin{1};  % index list or array
            end
         else
            range = obj.dataRangeFromVarargin(varargin{:});
            i = imfrmtRangeToIndex(obj.dataSize, obj.dataFormat, range);
         end
      end
      
      function i = cellIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            if nargin > 2
               i = imsub2ind(obj.cellSize, [varargin{:}]);
            else
               i = varargin{1};
            end
         else
            range = obj.cellRangeFromVarargin(varargin{:});
            i = imfrmtRangeToIndex(obj.cellSize, obj.cellFormat, range);
         end
      end
      
      function i = rawDataIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            if nargin > 2
               i = imsub2ind(obj.rawDataSize, [varargin{:}]);
            else
               i = varargin{1};
            end
         else
            range = obj.rawDataRangeFromVarargin(varargin{:});
            i = imfrmtRangeToIndex(obj.rawDataSize, obj.rawDataFormat, range);
         end
      end
      
      function i = rawCellIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            if nargin > 2
               i = imsub2ind(obj.rawCellSize, [varargin{:}]);
            else
               i = varargin{1};
            end
         else
            range = obj.rawCellRangeFromVarargin(varargin{:});
            i = imfrmtRangeToIndex(obj.rawCellSize, obj.rawCellFormat, range);
         end
      end
      
      
      function coords = dataIndexToCoordinate(obj, id, varargin)
         if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
            id = imsub2ind(obj.dataSize, [id, varargin{:}]);
         end 
         coords = imfrmtIndexToCoordinate(obj.dataSize, obj.dataFormat, id);
      end
      
      function coords = cellIndexToCoordinate(obj, id, varargin)
         if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
            id = imsub2ind(obj.cellSize, [id, varargin{:}]);
         end 
         coords = imfrmtIndexToCoordinate(obj.cellSize, obj.cellFormat, id);
      end
      
      
      function coords = rawDataIndexToCoordinate(obj, id, varargin)
         if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
            id = imsub2ind(obj.rawDataSize, [id, varargin{:}]);
         end 
         coords = imfrmtIndexToCoordinate(obj.rawDataSize, obj.rawDataFormat, id);
      end
      
      function coords = rawCellIndexToCoordinate(obj, id, varargin)
         if nargin >= 3 && isnumeric(varargin{1}) && isscalar(varargin{1})
            id = imsub2ind(obj.rawCellSize, [id, varargin{:}]);
         end 
         coords = imfrmtIndexToCoordinate(obj.rawCellSize, obj.rawCellFormat, id);
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
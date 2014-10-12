classdef ImageSource < ImageInfo
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use obj.data to return the actual image data
   %
   % required functions: datasize, dataformat, color, data
   %
   
   % Note: for translation to python, the cache structure should be via a global imge cache class
   %       that links ImageSource classes to the cached image
   %
   properties   
      icache = true;        % (optional) cache the data or not
      icelldata    = [];    % (cached) image data as cell array
      
      irawcache = true;     % (optional) cache the raw cell data
      irawcelldata = [];    % (cached) raw data
   end
        
   methods   
      function obj = ImageSource(varargin)  % basic constructor
         %
         % ImageSource()
         % ImageSource(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageSource') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj = obj.fromData(varargin{1});
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            obj = obj.fromParameter(varargin);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         disp source
         obj = classFromParameter(obj, 'i', varargin);
      end
      
      function obj = fromData(obj, data)
         fromData@ImageInfo(obj, data);
         if ~iscell(data)
            data = {data};
         end
         obj.irawcelldata = data; 
      end
   
   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % cache
      
      function cd = cellDataCache(obj)
         cd = obj.icelldata;
      end
      
      function obj = clearCellDataCache(obj)
         obj.icelldata = [];
      end
      
      function obj = setCellDataCaching(obj, c)
         if c ~= obj.icache
            obj.icache = c;
            obj.clearCellDataCache();
         end
      end

      function cd = rawCellDataCache(obj)
         cd = obj.irawcelldata;
      end
             
      function obj = clearRawCellDataCache(obj)
         obj.irawcelldata = [];
      end
      
      function obj = setRawCellDataCaching(obj, c)
         if c ~= obj.irawcache
            obj.irawcache = c;
            obj.clearRawCellDataCache();
         end
      end
      
      function obj = setCaching(obj,c)
         obj.setCellDataCaching(c);
         obj.setRawCellDataCaching(c);
      end
      
      
      function obj = cleatCache(obj)
         obj.clearCellDataCache;
         obj.clearRawCellDataCache;
      end
         
      
      
      function initializeCache(obj)
         % initialize cache by reading all images
         obj.setCaching(true);
         obj.icelldata = obj.cellData; % read all data in to cell
      end
      
      function initializeRawCache(obj)
         % initialize cache by reading all images
         obj.setCaching(true);
         obj.icelldata = obj.rawCellData; % read all data in to cell
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % data access methods

      function d = getDataFromRaw(obj, varargin)  % obtain the image data
         %
         % d = getData(obj, varargin)
         %
         % description:
         %     obtain data from raw data
         
         % range
         range = parseParameter(varargin);
         
         % check if all cell dims are singeltons
         cid = obj.cellIndex(range);
         if length(cid) > 1
            error('%s: getData: cell dimensions not specified to singeltons!', class(obj));
         end
         
         % find raw ranges necessary to load the raw data
         [rawRange, rawReshapeSize] = imfrmtReshapeInverseCellDataRange(obj.dataSize, obj.dataFormat, obj.rawDataFormat, ...
                                                                 obj.cellSize, obj.cellFormat, obj.rawCellFormat, obj.reshapeFrom, ...
                                                                 obj.reshapeTo, obj.reshapeSize, range);

         %rawRange
         %rawReshapeSize
                                                 
                                                              
         % load raw cell data
         cd = obj.getRawCellData(rawRange);
         
         % transform to output cell data 
         cd = imfrmtReshapeCellData(cd, obj.rawDataFormat, obj.rawCellFormat, obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, rawReshapeSize);
   
         % cd should be singeston
         if length(cd) ~= 1
            error('%s: getData: internal error, check get RawCellData routine.', class(obj));
         end
         
         d = cd{1};
      end

      function obj = setData(obj, d, varargin)
         id = obj.dataIndex(varargin{:});
         if length(id) ~= 1
            error('%s: setData: cell dimensions not specified to singeltons!', class(obj));
         end 
         obj.icelldata{id} = d;
      end
      
      
      function d = getRawData(obj, varargin)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         
         cid = obj.rawCellIndex(varargin{:});
         if length(cid) > 1
            error('%s: getRawData: cell dimensinos are not specified to be singeltons!', class(obj));
         end
         d = obj.irawcelldata{cid}; % trivial here
      end
      
      function obj = setRawData(obj, d, varargin)
         id = obj.rawDataIndex(varargin{:});
         if length(id) ~= 1
            error('%s: setRawData: cell dimensions not specified to singeltons!', class(obj));
         end 
         obj.irawcelldata{id} = d;
      end
      
      
      function d = getRawCellData(obj, varargin)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         
         if nargin > 1
            id = obj.rawCellIndex(varargin{:});
            d = obj.irawcelldata(id);     
         else
            d = obj.irawcelldata; % trivial here
         end
      end

      function obj = setRawCellData(obj, cd)  % set the image data
         obj.irawcelldata = cd;
      end

      function cd = getCellData(obj, varargin)
         cd = getCellDataFromRaw(obj, varargin{:});
      end
      
      function obj = setCellData(obj, cd)
         obj.icelldata = cd;
      end
      
      
      function cd = getCellDataFromRaw(obj, varargin)
         %
         % cd = getCellDataFromRaw(obj, varargin)
         %
         % description:
         %     obtain cell data from raw data and reshape specifications
         
         % range parameter
         range = imfrmtParseRangeLast(obj.irange, varargin);
         
         % find raw ranges necessary to load the raw data
         [rawRange, rawReshapeSize] = imfrmtReshapeInverseCellDataRange(obj.dataSize, obj.dataFormat, obj.rawDataFormat, obj.cellSize, obj.cellFormat, obj.rawCellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize, range);
         
         % load raw data as cell data
         cd = obj.getRawCellData(rawRange);
         
%          disp IS
%          var2char({ obj.rawDataFormat, obj.rawCellFormat, obj.dataFormat, obj.cellFormat})
%          size(cd)
%          size(cd{1})
         
         % transform to output cell data 
         cd = imfrmtReshapeCellData(cd, obj.rawDataFormat, obj.rawCellFormat, obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, rawReshapeSize);
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % data methods      

      function d = data(obj, varargin)
         if obj.icache % basic caching
            
            id = obj.cellIndex(varargin{:});
            if length(id) ~= 1
               error('%s: data: cell size not a singelton!');
            end
            
            if isempty(obj.icelldata)
               obj.icelldata = cell(obj.cellSize);
            end
               
            if ~isempty(obj.icelldata{id})
               d = obj.icelldata{id};
            else
               d = obj.getDataFromRaw(varargin{:});
               obj.icelldata{id} = d;
            end

         else
            d = obj.getDataFromRaw(varargin{:});
         end
      end
      
      
      function c = cell(obj, varargin)
         if obj.icache % basic caching
            
            id = obj.cellIndex(varargin{:});

            if isempty(obj.icelldata)
               obj.icelldata = cell(obj.cellSize);
            end

            idload = cellfun(@isempty, obj.icelldata(id));
            idload = id(idload);
               
            for i = 1:length(idload)
               obj.icelldata{idload(i)} = obj.data(idload(i));
            end
            
            c = obj.icelldata(id);
         else
            c = obj.getCellDataFromRaw(varargin{:});
         end
      end
      
          
      function d = dataSubset(obj, varargin)
         %
         % d = dataSubset(obj, datasepc)
         %
         % description:
         %    subset of the data given the data specifications datasepc
         %
         
         d = imfrmtDataSubset(obj.data, obj.dataFormat, varargin);
      end
      
      function d = dataExtract(obj, roi)
         %
         % d = subdata(obj, datasepc)
         %
         % description:
         %    extract a subset of the data given the spatial roi 
         %
         
         d = roi.extractData(obj.data); 
      end
      

       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% reshaping / reformatting data     
      
      function obj = setDataFormat(obj, dfrmt)         
         if ~isempty(obj.icelldata)
            obj.icelldata = obj.reformatData(obj.icelldata, obj.dataFormat, dfrmt);
         end
         setDataFormat@ImageInfo(obj, dfrmt);
      end
      
 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% casting and intensity scaling
     
      function m = initializeRawDataIntensity(obj, varargin)
         cd = obj.rawCellData(varargin);
         m = cellfun(@(x) min(x(:)), cd);
         obj.iminrawintensity = min(m(:));
         m = cellfun(@(x) max(x(:)), cd);
         obj.imaxrawintensity = max(m(:));  
      end
      
      function m = initializeDataIntensity(obj, varargin)
         cd = obj.cellData(varargin);
         m = cellfun(@(x) min(x(:)), cd);
         obj.imindataintensity = min(m(:));
         m = cellfun(@(x) max(x(:)), cd);
         obj.imaxdataintensity = max(m(:));   
      end 
      
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % information / visulaization 
      
      function plot(obj)
         plotImageSource(obj);
      end
      
   end
   
end
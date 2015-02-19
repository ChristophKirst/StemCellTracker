 classdef ImageSourceBF < ImageSource
   %
   % ImageSourceBF class provdes access to images in bio-format
   % 

   properties 
      ifilename       = '';  % bio-format file name
      ireaderInstance = [];  % instance of loci.formats.ChannelSeparator
   end
   
   properties (Dependent = true)
      ireader
   end

   methods
      
      function obj = ImageSourceBF(varargin) % constructor
         %
         % ImageSourceBF()
         % ImageSourceBF(filename)
         % ImageSourceBF(...,fieldname, fieldvalue,...)
         %
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceBF') %% copy constructor
               obj = copy(varargin{1});
            elseif ischar(varargin{1})
               obj.fromFile(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
         else
            obj.fromParmeter(varargin{:});
         end
         
         obj.icache = false;  % acces to images is almost as fast, also files are typically very big
         
         %bfinitialize;
      end
      
      function delete(obj)
         if ~isempty(obj.ireader)
            %fprintf('%s: closing file: %s\n', class(obj), obj.ifilename);
            obj.ireader.close();
         end
      end 
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, 'i', varargin{:});
      end
      
      function obj = fromFile(obj, filename, varargin)
         %
         % obj = fromFile(obj, filename, varargin)
         %
         % description:
         %      initialize source from a given file name

         if ~exist(filename, 'file')
            error('%s: fromFile: file %s does not exists!', class(obj), filename);
         end

         % initialize the reader
         obj.ifilename       = filename;
         
         % read from specific series if specified
         rawRange = obj.rawRange(varargin);
         if isfield(rawRange, 'S')
            s = {'S', rawRange.S};
         elseif isfield(rawRange, 's')
            s = {'s', rawRange.s};
         else
            s = {};
         end

         [info, obj.ireaderInstance] = imreadBFInfo(filename, s{:});
         obj.fromImageInfo(info);
      end
      
      function obj = initializeFromSeries(obj, s)
         obj.initializeFromFile(obj.ifilename, 'S', s);
      end
      
      function obj = fromImageInfo(obj, info)
         fromImageInfo@ImageInfo(obj, info);
      end

      function obj = initializeRangeFromFile(obj, varargin)
         obj.irange = obj.rawRangeFromFile(varargin{:});
      end

      function obj = initializeFromTiling(obj, varargin)
         % automatic tiling
         [ts, tf] = obj.tileSizeAndFormat(varargin{:});
         obj.setReshape('S', 'UV', ts);
         obj.setCellFormat(tf);
      end
      
      
      % for parallel processing the non-serializeable ChannerlReder is not broadcasted to he workers -> check if intialized
      function ir = get.ireader(obj)
         if isempty(obj.ireaderInstance)
            obj.ireaderInstance = imreadBFReader(obj.filename);
         end
         ir = obj.ireaderInstance;
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % bf functions
      
      function fn = currentFilename(obj)
         fn = char(obj.ireader.getCurrentFile());
      end
      
      function tp = type(obj)
         tp = char(obj.ireader.getFormat());
      end
      
      function fn = filename(obj)
         fn = obj.ifilename;
      end
      
      function ir = reader(obj)
         ir = obj.ireader;
      end
      
      function pos = stagePositions(obj, varargin)
         rgs = obj.rawRange(varargin{:});
         pos = imreadBFStagePositions(obj.ireader, rgs);
         %pos = imfrmtReshapeCellData(pos, 'XY', 'ZCTS', obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function vs = voxelSize(obj, varargin)
         vs = imreadBFVoxelSize(obj.ireader, varargin{:});
      end
      
      function ts = timeStamps(obj, varargin)
         ts = imreadBFTimeStamps(obj.ireader, varargin{:});
      end
      
      function et = exposureTimes(obj, varargin)
         et = imreadBFExposureTimes(obj.ireader, varargin{:});
      end
  
      function [tsi, tpos] = tileSize(obj,varargin)
         [tsi, tpos] = imreadBFTileSize(obj.ireader, varargin{:});
      end
      
      function [tsi, tfrmt] = tileSizeAndFormat(obj,varargin)
         [tsi, tfrmt] = imreadBFTileSizeAndFormat(obj.ireader, varargin{:});
      end

      function ts = waveLengths(obj,varargin)
         ts = imreadBFWaveLengths(obj.ireader, varargin{:});
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % ranges
      
      function tr = rawRangeFromFile(obj, varargin)
         tr = imreadBFRange(obj.ireader, varargin{:});
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access

      function d = getRawData(obj, range)
         % get raw data, varargin are range specs in raw format / size coordinates
    
         % check if all cell dims are singeltons
         cid = obj.rawCellIndex(range);
         if length(cid) > 1
            error('%s: getRawData: cell dimensions are not specified to be singeltons!', class(obj));
         end
        
         if obj.irawcache
            if isempty(obj.irawcelldata)
               cs = imfrmtAllocateSize(obj.rawCellSize);
               obj.irawceldata = cell(cs);
            end
            
            if isempty(obj.irawcelldata{cid}) 
               % read full data from file
               d = imreadBF(obj.ireader, obj.rawCellFormat, cid, 'squeeze', false, 'metadata', false);
               d = imfrmtReformat(d, 'XYZCT', obj.rawDataFormat); 
               
               % cach it
               obj.irawcelldata{cid} = d;
            else
               d = obj.irawcelldata{cid};
            end
               
            % reduce to range
            if ~isemptystruct(range)
               d = imfrmtDataSubset(d, obj.rawDataFormat, range);
            end
         else
            % read directly
            d = imreadBF(obj.ireader, range, 'squeeze', false, 'metadata', false);
            d = imfrmtReformat(d, 'XYZCT', obj.rawDataFormat);  
         end
      end
      
      
      function d = getRawCellData(obj, range) 
         % get raw cell data, varargin are range specs in raw format / size coordinates

         if obj.irawcache            
            if isempty(obj.irawcelldata)
               cs = imfrmtAllocateSize(obj.rawCellSize);
               obj.irawcelldata = cell(cs);
            end

            cid = obj.rawCellIndex(range);

            cidload = cellfun(@isempty, obj.irawcelldata(cid));
            cidload = cid(cidload);

            if ~isempty(cidload)
               % read missing data from file
               d = imreadBF(obj.ireader, obj.rawCellFormat, cidload , 'squeeze', false, 'metadata', false);
               if ~iscell(d)
                  d = {d};
               end
               d = imfrmtReformatCellData(d, 'XYZCT', obj.rawCellFormat, obj.rawDataFormat, obj.rawCellFormat);
               
               % cache
               obj.irawcelldata(cidload) = d;
            end
            
            % return requested cells
            d = obj.irawcelldata(cid);
            
         else
            d = imreadBF(obj.ireader, range, 'squeeze', false, 'metadata', false); 
            if ~iscell(d)
               d = {d};
            end
            d = imfrmtReformatCellData(d, 'XYZCT','S', obj.rawDataFormat, obj.rawCellFormat);
         end
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % plotting 
      
      function plot(obj, varargin)
         imgs = obj.cell(varargin{:});
         implot(imgs{1}, varargin{:});
      end
      
      function plotCells(obj, varargin)
         imgs = obj.cell(varargin{:});
         if nargin> 1 && isnumeric(varargin{1})
            varargin = varargin(2:end);
         end
         implottiling(imgs, varargin{:});
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % info
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj, 'BF');
         info = [info, '\nchannel:    '  var2char(obj.channelName)];
         info = [info, '\nfilename:   ', var2char(obj.filename)];
         info = [info, '\ntype:       ', var2char(obj.type)];
      end

   end
   
   
   methods(Access = protected)
      % Override copyElement method:
      function cpObj = copyElement(obj)
         %disp copyElement
         % Make a shallow copy of all four properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the DeepCp object
         cpObj.ireaderInstance = copy(obj.ireaderInstance);
      end
   end
   
end
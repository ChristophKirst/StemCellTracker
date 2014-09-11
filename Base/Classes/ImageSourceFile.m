classdef ImageSourceFile < ImageSource
   %
   % ImageSourceFile class represents image stored in a single file obtained by imread_bf(ifilename)
   % 

   properties 
      ifilename    = '';                % the file name of the image
      ireader      = @imread_bf;        % the reading routine used to 
      iinforeader  = @imread_bf_info    % the reading routine to obtain the info of the image
   end

   methods
      function obj = ImageSourceFile(varargin) % constructor
         %
         % ImageSourceFile()
         % ImageSourceFile(imagesourcefile)
         % ImageSourceFile(...,fieldname, fieldvalue,...)
         %

         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceFile') %% copy constructor
               obj = copy(varargin{1});
            elseif ischar(varargin{1})
               obj.ifilename = varargin{1};
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
         else
            for i = 1:2:nargin % constructor from arguments
               if ~ischar(varargin{i})
                  error('%s: invalid constructor input, expects char at position %g',class(obj), i);
               end
               if isprop(obj, lower(varargin{i}))
                  obj.(lower(varargin{i})) = varargin{i+1};
               else
                  warning('%s: unknown property name: %s ', class(obj), varargin{i})
               end
            end
         end
         
         if isempty(obj.ifilename)
            error('ImageSourceFile: filename needs to be specified');
         end
         
         if ~isfile(obj.ifilename)
            error('ImageSourceFile: filename %s does not point to a file', obj.ifilename);
         end

      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % access methods methods 
      % to obtain non-cached cached/preset info

      function img = getData(obj, varargin)
         img = obj.ireader(obj.ifilename);
      end
      
      function obj = setData(obj, varargin) % set the image data
         warning('ImageSourceFile: writing image data not suported, changes not propagated to disk!');
         obj.idata  = varargin{1}; 
         obj.chache = true;
         i = imdata2info(obj.idata);
         obj.setInfo(i);
      end

      function i = getInfo(obj, varargin)
         i = obj.iinforeader(obj.ifilename, varargin{:});
      end


      % infoString
      function istr = infoString(obj)
         istr = infoString@ImageSource(obj, 'File');
         istr = [istr, '\nfilename:   ', obj.ifilename];
         istr = [istr, '\nreader:     ', func2str(obj.ireader)];
         istr = [istr, '\ninforeader: ', func2str(obj.iinforeader)]; 
      end

   end
   
   
end
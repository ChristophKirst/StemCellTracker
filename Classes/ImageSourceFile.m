classdef ImageSourceFile < ImageSource
   %
   % ImageSourceFile class represents image stored a single file
   % 

   properties 
      filename = '';
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
               obj.filename = varargin{1};
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
                  warning('%s: unknown property name: %s ', class(obj), lower(varargin{i}))
               end
            end
         end
         
         if isempty(obj.filename)
            error('ImageSourceFile: filename needs to be specified');
         end
         
         if ~isfile(obj.filename)
            error('ImageSourceFile: filename does not point to a file');
         end
         
         obj.initialize()
         
      end
      
      
      function img = getRawData(obj, varargin)
         img = imread(obj.filename);
      end
      
      function info = readInfo(obj, varargin)
         info = imfinfo(obj.filename);
      end
      
      function info = getInfo(obj, varargin)
         basic = getInfo@ImageSource(obj);
         info = obj.readInfo();
         for f = fieldnames(basic)
            info.(f) = basic.(f);
         end
      end

      
      function siz = getSize(obj, varargin)
         info = obj.readInfo();
         
         z = length(info);
         w = info(1).Width;
         h = info(1).Height;
         
         if nargin < 2
            c = info(1).ColorType;
            switch c
               case 'grayscale'
                  c = 1;
               case 'truecolor'
                  c = 3;
               otherwise
                  warning('ImageSourceFile: cannot infer channel number, specify size directly');
                  c = 1;
            end
         else
            c = varargin{1};
         end
         
         if c == 1
            if z == 1
               siz = [w, h];
            else
               siz = [w,h,z];
            end
         else
            if z == 1
               siz = [w,h,c];
            else
               siz = [w,h,z,c];
            end
         end
         
      end
      
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj);
         info = [info, '\nfilename: ', obj.filename];
      end

   end
   
   
end
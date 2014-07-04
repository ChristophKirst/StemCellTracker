classdef ImageSourceReadCommand < ImageSource
   %
   % ImageSourcReadCommand class represents an image that is accessed by evaluating a command string
   % 

   properties 
      command = '';
   end

   methods
      function obj = ImageSourceReadCommand(varargin) % constructor
         %
         % ImageSourcReadCommand()
         % ImageSourcReadCommand(command)
         % ImageSourcReadCommand(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceReadCommand') %% copy constructor
               obj = copy(varargin{1});
            elseif ischar(varargin{1})
               obj.command = varargin{1};
            else
               error('ImageSourceReadCommand: invalid constructor input, expects char at position %g', 1);
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
         
         
         if isempty(obj.command)
            error('ImageSourceReadCommand: command needs to be specified in constructor!');
         end
         
         if isempty(obj.size)
            obj.size = obj.getSize();
         end

         if isempty(obj.format)
            obj.format = obj.getFormat();
         end

      end

      function img = getRawData(obj, varargin)
         img = eval(obj.command);
      end
      
      function siz = getSize(obj, varargin)        
         try
            res = eval(obj.command);
         catch
            warning('ImageSourceReadCommand: running command failed, cannot infer size');
            return
         end

         siz = size(res);
      end
      
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj);
         info = [info, '\ncommand: ', obj.command];
      end

   end
   
   
end
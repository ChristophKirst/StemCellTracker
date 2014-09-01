classdef Object < handle
%
% Object class storing basic information necessary for tracking objects
%
% See also: DataObject, Frame, Trajectory

   properties
      r = [0; 0];     % spatial position
      time = 0;       % time 
      
      volume    = 0;  % area or volume 
      intensity = 0;  % intensity 
      id        = 0;  % id
      
      type = [];      % object type ([] for none)
   end
   
   methods
      function obj = Object(varargin)  % simple constructor
      %
      % Object()
      % Object(...,fieldname, fieldvalue,...)
      %
         if nargin == 1 && isa(varargin{1}, 'Object') %% copy constructor 
            oldobj = varargin{1};
            nobjs = length(oldobj);
            obj(nobjs) = Object();
            for k =1:nobjs
               obj(k).r      = oldobj(k).r;
               obj(k).time   = oldobj(k).time;
               obj(k).volume = oldobj(k).volume;
               obj(k).id     = oldobj(k).id;
               obj(k).type   = oldobj(k).type;
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
      end

      
      function newobj = copy(obj)
         %
         % c = copy(obj)
         %
         % description:
         %     deep copy of the object
         %
         nobjs = length(obj);
         newobj(nobjs) = Object();
         for k =1:nobjs
            newobj(k).r      = obj(k).r;
            newobj(k).time   = obj(k).time;
            newobj(k).volume = obj(k).volume;
            newobj(k).id     = obj(k).id;
            newobj(k).type   = obj(k).type;
         end  
      end
      
  
      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % description:
         %     spatial dimension of the underlying coordinates
         %
            d = size(obj(1).r, 1);
      end
     
      
      function a = toArray(obj)
      %
      % data = toArray(obj)
      %
      % convert data to column vector
      %     
         a = [obj.time; obj.r; obj.volume; obj.intensity; obj.id; obj.type];      
      end
      
 
      function obj = transformCoordinates(obj, R, T, C)
      %
      % obj = transformCoordinates(obj, R, T, C)
      %
      % applies rotation R, scaling C  and translation T to coordinates r
      %  
         n = length(obj);
         trnsf = C * R * [ obj.r ] + repmat(T,1,n);
         for i=1:n
            obj(i).r = trnsf(:,i);
         end
      end

   end % methods
end % classdef
      
      
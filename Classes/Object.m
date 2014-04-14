classdef Object < handle
%
% Object class storing basic information necessary for tracking objects
%
% See also: Cell, Frame, Trajectory
   properties
      r = [0; 0];     % spatial position
      time = 0;       % time 
      
      volume    = 0;  % area or volume ([] for none)
      intensity = 0;  % intensity ([] for none)
      id        = 0;  % id ([] for none)
      
      type = [];      % object type ([] for none)
   end
   
   %properties (Dependent)
   %   length         % length of data entries obtained when using toArray
   %end
   
   methods
      function obj = Object(varargin)  % simple constructor
      %
      % Object()
      % Object(...,fieldname, fieldvalue,...)
      %
         for i = 1:2:nargin
            obj.(varargin{i}) = varargin{i+1};
         end
      end

      
      function c = copy(obj)
      %
      % c = copy(obj)
      %
      % deep copy of the object
      %   
         n = length(obj);
         for i=n:-1:1
            c(i) = Object(obj(i));
         end
         
      end
      
      function d = dim(obj)
         if ~isempty(obj)
            d = size(obj(1).r, 1);
         else
            error('Object: internal error!');
         end
      end
      
      
      %function l = get.length(obj)
      %      l = 1 + obj(1).dim + ~isempty(obj(1).volume) + ~isempty(obj(1).intensity) + ~isempty(obj(1).id) + size(obj(1).type,1);
      %end
      
      function data = toArray(obj)
      %
      % data = toArray(obj)
      %
      % convert data to column vector
      %    
         if length(obj) == 1
            data = [obj.time; obj.r; obj.volume;  obj.intensity; obj.id; obj.type];
         else
            
            d = [];
            if ~isempty(obj(1).volume) 
               d = [obj.volume];
            end
            if ~isempty(obj(1).intensity)
               if ~isempty(d)
                  d = vertcat(d, [obj.intensity]);
               else
                  d = [obj.intensity];
               end
            end

            if ~isempty(obj(1).id)
               if ~isempty(d)
                  d = vertcat(d, [obj.id]);
               else
                  d = [obj.id];
               end
            end 

            if ~isempty(obj(1).type)
               if ~isempty(d)
                  d = vertcat(d, [obj.type]);
               else
                  d = [obj.id];
               end
            end  

            if ~isempty(d)    
               data = vertcat( [obj.time],  [obj.r], d);
            else
               data = vertcat( [obj.time],  [obj.r]);
            end
         end       
      end % toArray
      
 
      function obj = transformCoordinates(obj, R, T, C)
      %
      % obj = transformCoordinates(obj, R, T, C)
      %
      % applies rotation R, scaling C  and translation T to coordinates r
      %  
         n = length(obj);
         
         if n==1
            obj.r = C * R * obj.r + T;
         else
            trnsf = C * R * [ obj.r ] + repmat(T,1,n);
            for i=1:n
               obj(i).r = trnsf(:,i);
            end
         end
         
      end
     
      

   end % methods
end % classdef
      
      
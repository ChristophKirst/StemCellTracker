function param = parseParameter(varargin)

if nargin  == 0
   param = [];
   return
end

if nargin >= 1
   if isstruct(varargin{1}) 
      param = varargin{1};
      vararg = varargin(2:end);
      narg = length(varargin)-1;
   else
      if nargin == 1
         if isempty(varargin{1})
            param = [];
            return
         else
            error('parseParameter: invalid parameter input, expects struct, [] or char at position %g', 1);
         end 
      else
         if isempty(varargin{1})
            param = struct();
            vararg = varargin(2:end);
            narg = length(varargin)-1;
         else
            vararg = varargin;
            narg = nargin;
            param = struct();
         end
      end
   end
end

for i = 1:2:narg % constructor from arguments
   if ~ischar(vararg{i})
      error('parseParameter: invalid constructor input, expects char at position %g', i);
   else
      param.(vararg{i}) = vararg{i+1};
   end
end

end
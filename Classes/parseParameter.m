function param = parseParameter(varargin)

if nargin == 0
   param = [];
   return
end

if nargin >= 1
   if isstruct(varargin{1}) 
      param = varargin{1};
      vararg = varargin(2:end);
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
         else
            vararg = varargin;
            param = struct();
         end
      end
   end
end

newparam = setParameter(vararg{:});
if ~isempty(newparam)
   parnames = fieldnames(newparam);

   for i = 1:length(parnames)
      param.(parnames{i}) = newparam.(parnames{i});
   end
end

end
function imgp = implotis(obj, varargin)
%
% img = implotis(obj, varargin);
%
%
% description:
%     takes information on image an plots it accordingly
%
% input:
%     obj     a ImageSource class or one of its derivatives
%
% output:
%    imgp     (optional) the composed image that is plotted


param = parseParameter(varargin);
%ov = getParameter(param, 'overlay', false);

% 
% if ~ov
%   cl = obj.datasizeC;
%  
%   for c = 1:cl
%      imsubplot(cl,1,c)
%      imcolormap(obj.color{c});
%      img = implot(obj.subdata('c', c));
%   end
%  
% else
   if obj.datasizeC == 1
      imcolormap(obj.color{1});
      img = implot(obj.data);
   elseif obj.datasizeC == 3 
      if isequal(obj.color, {'r', 'g', 'b'})
         img = implot(obj.data);
      else
         % TODO: handle special color settings here
         img = implot(obj.data);
      end
   else
      % TODO: handle non-standard channels here
      img = implot(obj.subdata('c', [1,2,3]));
   end
%end

if ~isempty(obj.name)
   set(gcf, 'Name', var2char(obj.name));
end

if nargout > 0
   imgp = img;
end

end
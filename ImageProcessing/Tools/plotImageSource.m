function imgp = plotImageSource(obj, varargin)
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


% TODO: make nice

if obj.dataSizeC == 1 && ~isempty(obj.color)
   imcolormap(obj.color{1});
   img = implot(obj.data);
elseif obj.dataSizeC == 3
   if isequal(obj.color, {'r', 'g', 'b'})
      img = implottiling(obj.data);
   else
      % TODO: handle special color settings here
      img = implot(obj.data);
   end
else
   % TODO: handle non-standard channels here
   img = implot(obj.dataSubset('C', [1,2,3]));
end

if ~isempty(obj.name)
   set(gcf, 'Name', var2char(obj.name));
end

if nargout > 0
   imgp = img;
end

end
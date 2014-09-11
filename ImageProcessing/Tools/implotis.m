function img = implotis(obj, varargin)
%

if obj.sizeC == 1
   imcolormap(obj.color{1});
   img = implot(obj.data);
elseif obj.sizeC == 3
   if isequal(obj.color, {'r', 'g', 'b'})
      img = implot(obj.data);
   else
      % TODO: handle special color settings here
      img = implot(obj.data);
   end
else
   % TODO: handle non-standard channels here
   img = implot(obj.data);
end

if isentry(obj.iinfo, 'iname')
   set(gcf, 'Name', var2char(obj.iinfo.iname));
end

end
function clearclass(cname)
%
%  clearclass(classname)
%
% description:
%    clears all instances of classes with class name cname

var = evalin('base', 'whos');
pos = strcmp({var.class}, cname);
names = {var(pos).name};

if ~isempty(names)
   cmd = 'clear(';
   for i = 1:length(names)-1
      cmd = [cmd, '''', names{i}, ''', ']; %#ok<AGROW>
   end
   cmd = [cmd, '''', names{end}, ''');'];
   evalin('base', cmd);
end

end




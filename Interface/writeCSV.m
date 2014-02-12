function success = writeCSV(filename, data, header)
%
% success = writeCSV(filename, data, header)
%
% writes data with header to csv file
%

fid = fopen(filename,'w');

if fid == -1
   success = -1;
   return
end

if nargin > 2
   if iscell(header)
      header = cell2mat(cellfun(@(x) [x, ', '], header, 'UniformOutput', false));
      header = header(1:end-2);
   end
   
   fprintf(fid,'%s\r\n',char(header));
end

fclose(fid);

dlmwrite(filename, data, '-append','delimiter',',','newline','pc');
   
success = 1;

end
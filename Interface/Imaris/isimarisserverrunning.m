function b = isimarisserverrunning()
%
% b = isimarisserverrunning()
%
% description:
%     tests if Imaris Server application is running
%
% output:
%     b    true if sever is running
%
% See also: imarissart

   b = 0;
   if ispc()
      [~, result] = system(...
         'tasklist /NH /FI "IMAGENAME eq ImarisServerIce.exe"');
      if strfind(result, 'ImarisServerIce.exe')
         b = 1;
         return;
      end

   elseif ismac()
      [~, result] = system('ps aux | grep ImarisServerIce');
      if strfind(result, this.mImarisServerExePath)
         b = 1;
         return;
      end

   else
      error('imarisstart: unsupported platform.');
   end
end
function ijplot5d(img, varargin)
%
% ijplot5d(img, param)
%
% description:
%   display 5d image in ImageJ, assumes image of format XYZCT
%
% input:
%   img     image(x,y,z,c,t)
%   param   (optional)  parameter struct with
%           .title      title of image
%           .class      array class in which data should be transferred ('uint8')
%           .write      write data to files and open in exteranl imagej for large data sets (false) 
%           .directory  directory to write files to in case write is true (tempdir)
%           .filename   filename prefix to write files to in case write is true (tempname)
%
% Notes: to develop further commands to feed to ImageJ turn on the recorder while doing 
%        things interactively...Plugins -> Macros -> Recorder
%
% modified from EDS 12/2014

param = parseParameter(varargin);
title = getParameter(param, 'title', []);
cls   = getParameter(param, 'class', 'uint8');
wrt   = getParameter(param, 'write', false);

isize = size(img);
if length(isize) > 5
   error('ijplot5d: inconsistent image dimensions: %s\n', var2char(isize));
end
isize = padright(isize, 5, 1);


if wrt % use file transfer
   
   tmpdir = getParameter(param, 'directory', fullfile(tempdir, 'ijplot5d'));
   mkdir(tmpdir)

   tmpfn  = getParameter(param, 'filename',  []);
   if isempty(tmpfn) 
      [~, tmpfn] = fileparts(tempname);
   end

   fnlst = cell(isize(3) * isize(5));
   
   img = imrecast(img, cls, 'rescale', 'full');
   
   % write the data
   k = 1;
   for t = 1:isize(5)
      for z = 1:isize(3)
         fn = fullfile(tmpdir, [tmpfn, sprintf('Z%.3dT%.4d.tif', z, t)]);
         imwrite(squeeze(img(:,:,z,:,t)), fn);
         fnlst{k} = fn;
         k = k + 1;
      end
   end
  
   % open in imagej externally
   lfn = fullfile(tmpdir, [tmpfn '_filelist.txt']);
   fid = fopen(lfn, 'w');
   for i = 1:length(fnlst)
      fprintf(fid, '%s\n', fnlst{i});
   end
   fclose(fid);

   mfn = fullfile(tmpdir, [tmpfn '_openmacro.txt']);
   fid = fopen(mfn, 'w');
   fprintf(fid, 'run("Stack From List...", "open=%s");', lfn);
   fprintf(fid, 'run("Stack to Hyperstack...", "order=xyczt(default) channels=%d slices=%d frames=%d display=Color");', 1, isize(3), isize(5));
   fclose(fid);

   ip = ijpath;
   
   cmd = sprintf('cd %s; java -jar ij.jar -macro "%s" &', ip, mfn);
   system(cmd)
   
   %iji = ij.IJ; 
   %iji.runMacro('Stack From List...', sprintf('open=%s', lfn))
   

   
   cln = getParameter(param, 'cleanup', true);
   if cln
      %fprintf('opening image in ImageJ, press key to continue and delete temporary files...\n')
      %pause
      waitfor(msgbox('opening image in ImageJ ! press key to continue and delete temporary files..', 'ijplot5d'))
      
      system(['rm -r ', tmpdir]);
   end

else    % use interface

   
   % transfrom image to QXY where Q = C Z T
   nz = size(img,3);
   nc = size(img,4);
   nt = size(img,5);
   
   img = imfrmtReformat(img, 'XYZCT', 'XYCZT');
   img = reshape(img, [size(img,1), size(img,2), nz*nc*nt]);
   img = imfrmtReformat(img, 'XYZ', 'ZXY');
   img = imrecast(img, cls, 'rescale', 'full');
   %size(img)
   
   if isempty(title)
      title = '';
   end
   imp = MImageJ.createImage(title, 'zxy', img); % fast creation of image data
   
   imp.show();
   imp.updateAndDraw();
   
   %if isize(3) == 3
   %    MIJ.run('Stack to RGB');
   %else % check for nz=1 or nt=1
   
   cmd = sprintf('order=xyczt(default) channels=%d slices=%d frames=%d display=Composite', nc, nz, nt);
   MImageJ.run('Stack to Hyperstack...', cmd);
   
   %end


end


end
 

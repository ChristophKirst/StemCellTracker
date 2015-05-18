function imwriteNRRD(filename, data, varargin)
%
% imwriteNRRD(filename, data, pixelspacing, origin, encoding)
%
% description:
%      writes 3d data array to nrrd file
%
% input:
%      filename     file name, e.g. 'filename.nrrd'
%      data         3D data array to export
%      param        parameter struct with entries
%                   .pixelspacing   voxel size ([1,1,1])
%                   .origin         origin of image ([0,0,0])
%                   .encoding       encoding, 'raw', 'ascii', 'gzip' ('raw')

% Based on nrrdWrite

param = parseParameter(varargin);

pixelspacing = getParameter(param, 'pixelspacing', [1,1,1]);
origin       = getParameter(param, 'origin', [0,0,0]);
encoding     = getParameterInNames(param, 'encoding', {'raw', 'ascii', 'gzip'}, 'raw');

[pathf, fname, ext] = fileparts(filename);
format=ext(2:end);

% TODO: check this
data = permute(data, [2 1 3]);

dims=(size(data));   
ndims=length(dims);


% The same with output format
format = lower(format);
if ~(isequal(format,'nhdr') || isequal(format,'nrrd')) 
   error('imwriteNRRD: expect file extension: nhdr or nrrd, found %s!', format);
end
    
% Header

fid = fopen(filename, 'wb');
fprintf(fid,'NRRD0004\n');      % NRRD type 4

% Type of variable we're storing in our file
mtype   = class(data);
outtype = setDatatype(mtype);

fprintf(fid,['type: ', outtype, '\n']);

fprintf(fid,['dimension: ', num2str(ndims), '\n']);

% orienation
if isequal(ndims, 2)
   fprintf(fid,'space: left-posterior\n');
elseif isequal (ndims, 3)
   fprintf(fid,'space: left-posterior-superior\n');
end

fprintf(fid,['sizes: ', num2str(dims), '\n']);

if isequal(ndims, 2)
   fprintf(fid,['space directions: (', num2str(pixelspacing(1)), ...
      ',0) (0,', num2str(pixelspacing(2)), ')\n']);
   fprintf(fid,'kinds: domain domain\n');
elseif isequal (ndims, 3)
   fprintf(fid,['space directions: (', num2str(pixelspacing(1)), ...
      ',0,0) (0,', num2str(pixelspacing(2)), ',0) (0,0,', ...
      num2str(pixelspacing(3)), ')\n']);
   fprintf(fid,'kinds: domain domain domain\n');
end

fprintf(fid,['encoding: ', encoding, '\n']);

[~,~,endian] = computer();

if (isequal(endian, 'B'))
   fprintf(fid,'endian: big\n');
else
   fprintf(fid,'endian: little\n');
end

if isequal(ndims, 2)
   fprintf(fid,['space origin: (', num2str(origin(1)),',', num2str(origin(2)),')\n']);
elseif isequal (ndims, 3)
   fprintf(fid,['space origin: (', num2str(origin(1)), ...
      ',',num2str(origin(2)),',', num2str(origin(3)),')\n']);
end

if (isequal(format, 'nhdr')) 
   fprintf(fid, ['data file: ', [fname, '.', encoding], '\n']);
   
   fclose(fid);
   if isequal(length(pathf),0)
      fid = fopen([fname, '.', encoding], 'wb');
   else
      fid = fopen([pathf, filesep, fname, '.', encoding], 'wb');
   end
else
   fprintf(fid,'\n');
end

writeData(fid, data, outtype, encoding);
fclose(fid);
end

% ========================================================================
% Determine the datatype --> From mtype (matlab) to outtype (NRRD) -->    
% ========================================================================
function datatype = setDatatype(metaType)

% Determine the datatype
switch (metaType)
 case {'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64',...
       'uint64', 'double'}
   datatype = metaType;
  
 case {'single'}
  datatype = 'float';
  
 otherwise
  assert(false, 'Unknown datatype')
end

end
   

% ========================================================================
% writeData -->
% fidIn is the open file we're overwriting
% matrix - data that have to be written
% datatype - type of data: int8, string, double...
% encoding - raw, gzip, ascii
% ========================================================================
function ok = writeData(fidIn, matrix, datatype, encoding)

switch (encoding)
 case {'raw'}
  
  ok = fwrite(fidIn, matrix(:), datatype);
  
 case {'gzip'}
     
     % Store in a raw file before compressing
     tmpBase = tempname(pwd);
     tmpFile = [tmpBase '.gz'];
     fidTmpRaw = fopen(tmpBase, 'wb');
     assert(fidTmpRaw > 3, 'Could not open temporary file for GZIP compression');
     
     fwrite(fidTmpRaw, matrix(:), datatype);
     fclose(fidTmpRaw);
     
     % Now we gzip our raw file
     gzip(tmpBase);
     
     % Finally, we put this info into our nrrd file (fidIn)
     fidTmpRaw = fopen(tmpFile, 'rb');
     tmp = fread(fidTmpRaw, inf, [datatype '=>' datatype]);
     cleaner = onCleanup(@() fclose(fidTmpRaw));
     ok = fwrite (fidIn, tmp, datatype);
     
     delete (tmpBase);
     delete (tmpFile);


 case {'ascii'}
  
  ok = fprintf(fidIn,'%u ',matrix(:));
  %ok = fprintf(fidIn,matrix(:), class(matrix));
  
 otherwise
  assert(false, 'Unsupported encoding')
end

end


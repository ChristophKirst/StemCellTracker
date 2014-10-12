function imwrite_tiff(data, filename, varargin)
%
% imwrite_tiff(data, filename, varargin)
%
% description:
%    writes data to tiff file
% 
% input:
%    data      image data
%    filename  image filename
%    varargin  (optional) pair tagname, tagvalue

% See also: Tiff

t = Tiff(filename, 'w');

% Gray, RGB, GrayAlpha or RGBA

siz = size(data);
nc = size(data,3);

switch nc
   case 1
      t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
   case 2
      t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
      t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
   case 3
      t.setTag('Photometric', Tiff.Photometric.RGB);
   case 4
      t.setTag('Photometric', Tiff.Photometric.RGB);
      t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
   otherwise
      error('imwrite_tiff: image data is not Gray, RGB, GrayAlpha or RGBA');
end

% check data type

cls = class(data);

switch cls
   case 'uint8'
      t.setTag('BitsPerSample',8);
      t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
   case 'uint16'
      t.setTag('BitsPerSample',16);
      t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
   case 'uint32'
      t.setTag('BitsPerSample',32);
      t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
   case 'uint64'
      t.setTag('BitsPerSample',64);
      t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
   case 'int8'
      t.setTag('BitsPerSample',8);
      t.setTag('SampleFormat',Tiff.SampleFormat.Int);
   case 'int16'
      t.setTag('BitsPerSample',16);
      t.setTag('SampleFormat',Tiff.SampleFormat.Int);
   case 'int32'
      t.setTag('BitsPerSample',32);
      t.setTag('SampleFormat',Tiff.SampleFormat.Int);
   case 'int64'
      t.setTag('BitsPerSample',64);
      t.setTag('SampleFormat',Tiff.SampleFormat.Int);
   case {'single', 'double'}
      % rescale to uint16
      warning('imwrite_tiff: rescaling double image data to unit16');
      data = imrescale(data, 'class', 'uint16');
      t.setTag('BitsPerSample',16);
      t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
   otherwise
      error('imwrite_tiff: unknonw image data class %s for export!', cls);
end

t.setTag('ImageLength',siz(1));
t.setTag('ImageWidth', siz(2));
t.setTag('SamplesPerPixel',nc);

% general tags 
t.setTag('Compression',Tiff.Compression.None);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);

% set additional tags

tifftags = parseParameter(varargin{:});
if ~isempty(tifftags)
   for f = fieldnames(tifftags)'
      t.setTag(f{1}, tifftags.(f{1}));
   end
end

% write data
t.write(data);
t.close();

end


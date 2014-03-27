function success = writePAL(filename, cmap)
%
% success = writePAL(filename, cmap)
%
% description:
%     writes the rgb color map cmap to a pallete file pal (for import in
%     other software such as Imaris)
%
% input:
%     filename    name of pal file to save color map 
%     cmap        colormap 
%  
% output:
%     success     true if writing was succesful
%
% See also: colormap

 
%file specs
mf='n';							% Machine format
depth=size(cmap,1);			% Colormap depth
hlen=24;							% Header length
flen=hlen+(4*depth);
	
% open file
fid=fopen(filename,'w',mf);
if(fid<0)
   success = -1;
   return
end

% write RIFF signature
fwrite(fid,'RIFF','uint8',0,mf);
	
% write file length
fwrite(fid,flen-8,'uint32',0,mf);                        % 8 byte header (RIFF header)
	
% write PAL signature
fwrite(fid,'PAL ','uint8',0,mf);
	
% write data signature
fwrite(fid,'data','uint8',0,mf);

% write data block size
fwrite(fid,flen-20,'uint32',0,mf);                       % 20 byte header (RIFF + Chunk)

% write version number
fwrite(fid,[0,3],'uint8',0,mf);                          % Always 3

% write palette length
fwrite(fid,depth,'uint16',0,mf);

% write palette data
fwrite(fid,[cmap.*255,zeros(depth,1)]','uint8',0,mf);	   % RGBA tuples

% close file
fclose(fid);

end

   
   
   
   	
% 	% GUI path selection (if no path argument or directory)
% 	exts={'*.pal','Binary Palette File (*.pal)';'*.*','All Files (*.*)'};
% 	if(isempty(varargin))
% 		[file_name,file_dir,filter_index]=uiputfile(exts,'Save palette file');					%% No argument
% 		path=[file_dir,file_name];
% 	elseif(isdir(varargin{1}))
% 		[file_name,file_dir,filter_index]=uiputfile(exts,'Save palette file',varargin{1});		%% Directory argument
% 		path=[file_dir,file_name];
% 
%  % Validate manual path extension
% 	else
% 		path=varargin{1};																		%% Path argument
% 		[file_dir,file_name,file_ext]=fileparts(path);		
% 		switch(file_ext)
% 			case ''
% 				filter_index=2;
% 			case '.pal'
% 				filter_index=1;
% 			otherwise
% 				warning('cmap2pal:FileExt','File extension not .pal (%s). Changing extension to .pal',file_ext);
% 				path=[file_dir,file_name];
% 				filter_index=2;		
% 		end
% 		
% 	end
% 	
% % Catch cancel from uiputfile
% 	if(filter_index==0)
% 		return

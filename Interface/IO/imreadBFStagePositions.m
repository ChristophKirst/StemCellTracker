function pos = imreadBFStagePositions(name, varargin)
%
% pos = imreadBFStagePositions(name, varargin)
%
% decription:
%    tries to infer stage positions from bio-format image
%
% input:
%    name    filename or reader
%    param   parameter struct with entries
%            .S / .s  series ids
%            .T / .t  time ids
%            .C / .c  channel ids
%            .Z / .z  depth ids 
%             
% output:
%    pos     cell array of stage positions for each series
%            each series position is a cell with positions at z,c,t,s
%            positions are (x,y(,z) as array)
%
% note: 
%    its assumed that if multiple series are given the data sizes are the same as the first series 

param = parseParameter(varargin);

ireader = imreadBFReader(name);

[sids, nSeries] = imfrmtParseRangeToIndex(ireader.getSeriesCount(), 'S', param);
  
metadataStore = ireader.getMetadataStore();

se = ireader.getSeries();

ireader.setSeries(sids(1)-1);
[~, nTimes]  = imfrmtParseRangeToIndex(ireader.getSizeT(), 'T', param);
[~, nChls]   = imfrmtParseRangeToIndex(ireader.getSizeC(), 'C', param);
[~, nZs]     = imfrmtParseRangeToIndex(ireader.getSizeZ(), 'Z', param);

pos = cell(nZs, nChls, nTimes, nSeries);

sc = 1;
for s = sids - 1

   ireader.setSeries(s);
      
   tids  = imfrmtParseRangeToIndex(ireader.getSizeT(), 'T', param);
   cids  = imfrmtParseRangeToIndex(ireader.getSizeC(), 'C', param);
   zids  = imfrmtParseRangeToIndex(ireader.getSizeZ(), 'Z', param);

   zids = zids - 1;
   tids = tids - 1;
   cids = cids - 1;

   nPlanes = metadataStore.getPlaneCount(s);

   for p = 0:(nPlanes-1)
      
      z = find(zids == metadataStore.getPlaneTheZ(s,p).getValue, 1, 'first');
      c = find(cids == metadataStore.getPlaneTheC(s,p).getValue, 1, 'first');
      t = find(tids == metadataStore.getPlaneTheT(s,p).getValue, 1, 'first');
      
      if ~isempty(z) && ~isempty(t) && ~isempty(c)

         % Get X and Y position from current series metadata
         try
            xpos = metadataStore.getPlanePositionX(s,p).doubleValue();
         catch
            xpos = [];
         end

         try
            ypos = metadataStore.getPlanePositionY(s,p).doubleValue();
         catch
            ypos = [];
         end
         
         try
            zpos = metadataStore.getPlanePositionZ(s,p).doubleValue();
         catch
            zpos = [];
         end
      
         pos{z,c,t, sc} = [xpos, ypos, zpos]';
         
      end
   end
        
   sc = sc + 1;
       
end

% restore actual series
ireader.setSeries(se);

end


function cns = imreadBFChannelNames(name, varargin)
%
%  cns = imreadBFChannelNames(ireader, varargin)
%
% decription:
%    tries to infer channel names from bio-format image
%
% input:
%    name    filename or reader
%    param   parameter struct with entries
%            .S / .s  series ids
%             
% output:
%    cns     channel names
%

param = parseParameter(varargin);

ireader = imreadBFReader(name);

[sids, nSeries] = imfrmtParseRangeToIndex(ireader.getSeriesCount(), 'S', param);

metadataStore = ireader.getMetadataStore();

se = ireader.getSeries();

ireader.setSeries(sids(1)-1);
[cids, nColors] = imfrmtParseRangeToIndex(ireader.getSizeC(), 'C', param);

cns = cell(nColors,nSeries);

sc = 1;
for s = sids
          
   cc = 1;
   for ch = cids
 
      %try
         cn =char(metadataStore.getChannelName(s - 1, ch - 1));
      %catch
      %   name = [];
      %end
      if isempty(name)
         name = '';
      end
      
      cns{cc, sc} = cn;
      cc = cc + 1;
   end
   
   sc = sc + 1;
end
    
% restore actual series
ireader.setSeries(se);
  
end
   


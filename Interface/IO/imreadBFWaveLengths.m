function [emwl, exwl] = imreadBFWaveLengths(name, varargin)
%
% [emwl, exwl] = imreadBFWaveLengths(name, param)
%
% decription:
%    tries to infer excitation and emission wavelengths from bio-format image
%
% input:
%    name    filename or reader
%    param   parameter struct with entries
%            .n / .series  series ids
%            .c / .channel channel ids
%             
% output:
%    emwl    emission wave length
%    exwl    excitation wave length
%

param = parseParameter(varargin);

ireader = imreadBFReader(name);

[sids, nSeries] = imfrmtParseRangeToIndex(ireader.getSeriesCount(), 'S', param);

metadataStore = ireader.getMetadataStore();

se = ireader.getSeries();

ireader.setSeries(sids(1)-1);
[~, nColors] = imfrmtParseRangeToIndex(ireader.getSizeC(), 'C', param);

emwl = cell(nColors, nSeries);
exwl = cell(nColors, nSeries);

sc = 1;
for s = sids

   cc = 1;
   for ch = cids
      
      % Emission wavelength
      try
         em = metadataStore.getChannelEmissionWavelength(s - 1, ch - 1);
      catch
         em = [];
      end
      if isempty(em)
         em = NaN;
      else
         em = em.getValue();
      end
      emwl(cc, sc) = em;
      
      % Excitation wavelength
      try
         exc = metadataStore.getChannelExcitationWavelength(s - 1, ch - 1);
      catch
         exc = [];
      end
      if isempty(exc)
         exc = NaN;
      else
         exc = exc.getValue();
      end
      exwl(cc, sc) = exc;
      
      cc = cc + 1;
   end

   sc = sc + 1;
   
end

    
% restore actual series
ireader.setSeries(se);


end
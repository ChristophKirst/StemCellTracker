function outRange = imfrmtReshapeRangeOnce(inSize, inFrmt, outFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
% helper for imfrmtReshapeRange

% does the following:
%    ranges not being reshaped are left untouched
%    reshapes the ranges and checks if this gives a new valid multiplicative range
%    extra dimensions which are lost will be set to the first value in range if exists

% reshaping is trivial for empty range
if isemptystruct(range)
   outRange = range;
   return
end

% some checks
if length(reshapeSize) ~=length(reshapeTo)
   error('imfrmtReshapeRange: reshape format %s inconsistent with reshape size %s', reshapeTo, var2char(reshapeSize));
end

% range names
rangeNames = fieldnames(range);
rangeNamesl = lower(rangeNames);

posRangeFrom = ismember(rangeNamesl, num2cell(lower(reshapeFrom)));

outRange = range;

 % reshape
if any(posRangeFrom)
 
   % ranges not involved in reshaping
   range = rmfield(range, rangeNames(~posRangeFrom));
   %rangeNames = rangeNames(posRangeFrom);
   rangeNamesl = rangeNamesl(posRangeFrom);
   
   % supplement range with full missing reshape formats
   posFromRange = ismember(num2cell(lower(reshapeFrom)), rangeNamesl);
   sizeFrom = imfrmtReformatSize(inSize, inFrmt, reshapeFrom);

   ids = find(~posFromRange);
   for i = 1:length(ids)
      range.(reshapeFrom(ids(i))) = 1:sizeFrom(ids(i));
   end
   
   range = imfrmtReformatRange(sizeFrom, reshapeFrom, reshapeFrom, range);
   
   % indices
   idx = imfrmtRangeToIndex(sizeFrom, reshapeFrom, range);
   
   % to reshapeTo range
   %reshapeSize
   %reshapeTo
   %idx
   range = imfrmtRangeFromIndex(reshapeSize, reshapeTo, idx);
   
   
   % reformat to outFrmt
   posToRange = imfrmtPosition(outFrmt, reshapeTo);
   if length(posToRange) ~=length(reshapeTo)
      error('imfrmtReshapeRange: output format %s inconsistent with reshape format %s', outFrmt, reshapeTo);
   end

   range = imfrmtReformatRange(reshapeSize, reshapeTo, outFrmt(posToRange), range);

   % remove reshapeFrom and add reshapeTo ranges
   outRange = imfrmtRemoveRange(outRange, reshapeFrom);
   outRange = imfrmtRangeFromVarargin(outRange, range);
end

% check if dimensions are lost and if so set range to first 
fullOutFrmt = imfrmtReshapeFormat(inFrmt, reshapeFrom, reshapeTo);

[extrafrmt,extrapos] = imfrmtExtraFormats(outFrmt, fullOutFrmt);

if ~isempty(extrafrmt)
   % size of extra dims
   
   extrasize = imfrmtReshapeSize(inSize, inFrmt, fullOutFrmt, reshapeFrom, reshapeTo, reshapeSize);
   extrasize = extrasize(extrapos);
   
   if any(extrasize > 1)
      warning('imfrmtReshapeRange: input dimensions %s of size %s lost and set to singelton!', extrafrmt, var2char(extrasize))

      rgs = [num2cell(extrafrmt); num2cell(ones(1,length(extrafrmt)))];
      rgs = imfrmtRangeFromVarargin(struct(rgs{:}), outRange);
      rgs = imfrmtRangeFromSizeFormatAndVarargin(extrasize, extrafrmt, rgs);
      rgsn = fieldnames(rgs);
      for i = 1:length(rgsn)
         v = rgs.(rgsn{i});
         rgs.(rgsn{i}) = v(1);
      end
      outRange = imfrmtRangeFromVarargin(outRange, rgs);
   end
end




end

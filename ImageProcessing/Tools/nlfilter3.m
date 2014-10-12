function b = nlfilter3(varargin)
%NLFILTER3 General sliding-neighborhood operations.
%
% This is a simple extension of the nlfilter function from the Image
% Processing Toolbox to 3D matrix stacks.
%
%   B = NLFILTER3(A,[M N P],FUN) applies the function FUN to each
%   M-by-N-by-P sliding block of A.  FUN is a function that accepts an
%   M-by-N-by-P matrix as input and returns a scalar:
%
%       C = FUN(X)
%
%   FUN must be a FUNCTION_HANDLE.
%
%   C is the output value for the center pixel in the M-by-N/by-P block
%   X. NLFILTER3 calls FUN for each pixel in A. NLFILTER3 zero pads the
%   M-by-N-by-P block at the edges, if necessary.
%
%   B = NLFILTER3(A,'indexed',...) processes A as an indexed image, padding
%   with ones if A is of class double and zeros if A is of class uint8.
%
%   Class Support
%   -------------
%   The input image A can be of any class supported by FUN. The class of B
%   depends on the class of the output from FUN.
%
%   Aaron Ponti, 2008/08/11
%
%   NLFILTER: Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 5.20.4.6 $  $Date: 2006/06/15 20:09:14 $

[a, nhood, fun, params, padval] = parse_inputs(varargin{:});

% Expand A
[ma,na,pa] = size(a);
aa = mkconstarray(class(a), mean( a( : ) ), size(a)+nhood-1);
aa( ...
    floor((nhood(1)-1)/2)+(1:ma), ...
    floor((nhood(2)-1)/2)+(1:na), ...
    floor((nhood(3)-1)/2)+(1:pa) ) = a;

% Find out what output type to make.
rows   = 0:(nhood(1)-1); 
cols   = 0:(nhood(2)-1);
planes = 0:(nhood(3)-1);
b = mkconstarray(class(feval(fun,aa(1+rows,1+cols,1+planes),params{:})), 0, size(a));
    
% Apply m-file to each neighborhood of a
f = waitbar(0,'Applying neighborhood operation...');
nOps = ma * na * pa; counter = 0;
for i=1:ma,
  for j=1:na,
      for k=1:pa,
          x = aa(i+rows,j+cols,k+planes);
          b(i,j,k) = feval(fun,x,params{:});
          counter = counter + 1;
          waitbar( counter / nOps );
      end
  end
end
close(f)

%%%
%%% Function parse_inputs
%%%
function [a, nhood, fun, params, padval] = parse_inputs(varargin)

switch nargin
case {0,1,2}
    eid = sprintf('Images:%s:tooFewInputs',mfilename);
    msg = 'Too few inputs to NLFILTER3';
    error(eid,'%s',msg); 
case 3
    if (strcmp(varargin{2},'indexed'))
        eid = sprintf('Images:%s:tooFewInputsIfIndexedImage',mfilename);
        msg = 'Too few inputs to NLFILTER3';
        error(eid,'%s',msg);
    else
        % NLFILTER3(A, [M N P], 'fun')
        a = varargin{1};
        nhood = varargin{2};
        fun = varargin{3};
        params = cell(0,0);
        padval = 0;
    end
    
otherwise
    if (strcmp(varargin{2},'indexed'))
        % NLFILTER3(A, 'indexed', [M N P], 'fun', P1, ...)
        a = varargin{1};
        nhood = varargin{3};
        fun = varargin{4};
        params = varargin(5:end);
        padval = 1;
        
    else
        % NLFILTER3(A, [M N P], 'fun', P1, ...)
        a = varargin{1};
        nhood = varargin{2};
        fun = varargin{3};
        params = varargin(4:end);
        padval = 0;
    end
end

if (isa(a,'uint8') || isa(a,'uint16'))
    padval = 0;
end

fun = fcnchk(fun,length(params));

function out = mkconstarray(class, value, size)
out = repmat(feval(class, value), size);

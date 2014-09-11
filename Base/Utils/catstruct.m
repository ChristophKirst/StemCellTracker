function s = catstruct(varargin)
%
%  s = catstruct(s1, s2, ...)
%  s = catstruct(s1, s2, ..., 'sorted')
%
% description:
%    concatenates two struct, overwriting values from the first
%    with the last but keeping the orderin of the first if possible
%
% input:
%    s1, s2,...   input structs to be concatenated
%    'sorted'     sort entries alphabetically
%    
% output:
%    s            concatenated struct
% 

% Note: based on catstruct by XXX??


narginchk(1,Inf) ;
N = nargin ;

if ~isstruct(varargin{end}),
    if isequal(varargin{end},'sorted'),
        sorted = 1 ;
        N = N-1 ;
        if N < 1
           narginchk(1,nargin-1) % generate error
        end
    else
        error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted".') ;
    end
else
    sorted = 0 ;
end

sz0 = [] ; % used to check that all inputs have the same size

% used to check for a few trivial cases
NonEmptyInputs = false(N,1) ; 
NonEmptyInputsN = 0 ;

% used to collect the fieldnames and the inputs
FN = cell(N,1) ;
VAL = cell(N,1) ;

% parse the inputs
for ii=1:N,
    X = varargin{ii} ;
    if ~isstruct(X),
        error('catstruct:InvalidArgument',['Argument #' num2str(ii) ' is not a structure.']) ;
    end
    
    if ~isempty(X),
        % empty structs are ignored
        if ii > 1 && ~isempty(sz0)
            if ~isequal(size(X), sz0)
                error('catstruct:UnequalSizes','All structures should have the same size.') ;
            end
        else
            sz0 = size(X) ;
        end
        NonEmptyInputsN = NonEmptyInputsN + 1 ;
        NonEmptyInputs(ii) = true ;
        FN{ii} = fieldnames(X) ;
        VAL{ii} = struct2cell(X) ;
    end
end

if NonEmptyInputsN == 0
    % all structures were empty
    s = struct([]) ;
elseif NonEmptyInputsN == 1,
    % there was only one non-empty structure
    s = varargin{NonEmptyInputs} ;
    if sorted,
        s = orderfields(s) ;
    end
else
    % there is actually something to concatenate
    FN = cat(1,FN{:}) ;    
    VAL = cat(1,VAL{:}) ;    
    FN = squeeze(FN) ;
    VAL = squeeze(VAL) ;
    
    
    %[UFN,ind] = unique(FN, 'last') ;
    if sorted
      [UFN,ind] = unique(FN, 'first') ; 
    else
      [UFN,ind] = unique(FN, 'stable') ;
    end

    if numel(UFN) ~= numel(FN),
%        warning('catstruct:DuplicatesFound','Fieldnames are not unique between structures.') ;
        sorted = 1 ;
    end
    
    if sorted,
        VAL = VAL(ind,:) ;
        FN = FN(ind,:) ;
    end
    
    s = cell2struct(VAL, FN);
    s = reshape(s, sz0) ; % reshape into original format
end

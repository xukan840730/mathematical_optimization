function options = emdoptimset(varargin)


%   EMDOPTIMSET returns a listing of the fields in the options structure as
%   well as valid parameters and the default parameter. It is written
%   according to the GAOPTIMSET (genetic algoritms toolbox) and OPTIMSET (optimization toolbox) matlab functions.
%   
%   OPTIONS = EMDOPTIMSET(@DI,'PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a structure with the
%   default parameters used in the Doubly-Iterative EMD (see ref. [1]) for all PARAM not specified, and will use the
%   passed arguments VALUE1, VALUE2, etc. for the specified PARAM1, PARAM2,
%   etc.
%
%   OPTIONS = EMDOPTIMSET(@standard,'PARAM1',VALUE1,'PARAM2',VALUE2,....) or just
%   OPTIONS = EMDOPTIMSET('PARAM1',VALUE1,'PARAM2',VALUE2,....) creates a structure with the
%   default parameters used in the standard EMD (see ref. [2]) for all PARAM not specified, and will use the
%   passed arguments VALUE1, VALUE2, etc. for the specified PARAM1, PARAM2, etc.
%
%   OPTIONS = EMDOPTIMSET(OLDOPTS,'PARAM1',VALUE1,'PARAM2',VALUE2,...) will
%   reassign those fields in OLDOPTS specified by PARAM1, PARAM2, ... to 
%   VALUE1, VALUE2, ...
%
%EMDOPTIMSET PARAMETERS
%
% Stopping:    - Several Stopping Criteria.
%              'Nsifts': Each IMF is computed with N siftings.
%              'Nsifts_AND_IMF': Each IMF is computed with at least N
%                   siftings. If after N siftings the IMF criterion is not
%                   fulfilled, i.e. the number of extrema is not equal or +1 to
%                   the number of zerrocrossings, then the siftings continue
%                   until the IMF criterion is become true.
%               'Nsifts_after_IMF': The IMF is computed using N extra
%                   siftings after the IMF criterion become true.
%               'Single_T': Single threshold stopping criterion (see ref.
%               [2]).
%               'Double_T': Double threshold stopping criterion (see ref.
%               [3]).
%
% Nodes:        - The method used for interpolation points selection
%               'Extrema': The maxima and minima used in the conventional
%                   EMD (i.e. as in ref. [2]).
%               'DI_Extrema': The interpolation points are computed
%                   according to the Doubly-iterative EMD principle (see ref.
%                   [1]).
%
% N:            - Positive integer used as 'N' in the Stopping parameter
% 
% T:            - Parameter 'T' in the Stopping parameter 'Single_T'
%                       (positive scalar, i.e. T=0.0001) and the parameter
%                       'Double_T' (three values vector, i.e.
%                       T=[0.05,0.5,0.05] (see ref [3]). 
%
% Order:        - Spline interpolation order (parameter 'q' in ref. [1] page 7). It should take odd integer numbers. (The Matlab Spline toolbox is
%                   needed for interpolation order higher than 3).
%
% in_N:         - Number of internal sifting iterations (parameter 'it' in ref. [1] page 7. It is used with
%                   the DI_Extrema only).
%
% in_Nskip:     - Number of external sifting iterations that the internal
%                   sifting iterations are skipped (parameter 'ex' in ref.
%                   [1] page 7. It is used with the DI_Extrema only).
%
% in_Order:     - Internal sifting Spline interpolation order (parameter
%                   'q' in ref. [1] page 7). It should take odd integer
%                   numbers. (The Matlab Spline toolbox is needed for
%                   interpolation order higher than 3). 
%
% Differentiation: - Method for estimating the derivatives.
%                    'diff_2p': For 2 points estimation.
%                    'diff_5p': For 5 points estimation.
%
% MaxN:         - Maximum Number of sifting iterations independently on
%                   whether the stopping criteria have been fulfilled or not
%
% IMFs:         - Maximum Number of extracted IMFs
%
% Disp:         - 1 for displaying the progress of the decomposition. 0 for
%                   "silent mode".
%
% EXAMPLES:
% options = emdoptimset(@di,'Stopping','Double_T','T',[0.05,0.5,0.05])
% options = emdoptimset(options,'IMFs', 1)
% The Parameters are not case sensitive, and it is not necessary to
% indicated with their full name. i.e.,
% options = emdoptimset('Differentiation','diff_2p','MaxN',100) is equivalent to
% options = emdoptimset('dif','diff_2p','Max',100)
%
% References
% [1] Yannis Kopsinis and Steve McLaughlin, “Improved EMD Using
% Doubly-Iterative Sifting and High Order Spline Interpolation,” EURASIP
% Journal on Advances in Signal Processing, vol. 2008, Article ID 128293.
%
% [2] N. E. Huang, Z. Shen, S. R. Long, et al., “The empirical mode
% decomposition and the Hilbert spectrum for nonlinear and non-stationary
% time series analysis,” Proceedings of the Royal Society of London A, vol.
% 454, no. 1971, 903 pages, 1998.   
%
% [3] G. Rilling, P. Flandrin, and P. Gonçalvès, “On empirical mode
% decomposition and its algorithms,” in Proceedings of the 6th IEEE/EURASIP
% Workshop on Nonlinear Signal and Image Processing (NSIP '03), Grado,
% Italy, June 2003.   
%
% Ver. 1.0
% Yannis Kopsinis, 23/04/2008
% kopsinis@ieee.org

optionsstandart=struct('Stopping','Nsifts_after_IMF',...
            'N',10,...
            'T',[],...
            'MaxN',5000,...
            'IMFs',50,...
            'Disp',1,...
            'Order',3,...
            'Nodes','Extrema',...
            'Differentiation','diff_5p',...
            'in_N',[],...
            'in_Nskip',[],...
            'in_Order',[]);

Names = fieldnames(optionsstandart);

names = lower(Names);
      
        
if nargin==0 && nargout==1
    options=optionsstandart;
    return
elseif nargin==0 && nargout==0
    fprintf('       Stopping: [''Nsifts''  |  ''Nsifts_AND_IMF''  |  {''Nsifts_after_IMF''}  |  ''Single_T''  |  ''Double_T''] \n')
    fprintf('          Nodes: [{''Extrema''}  |  ''DI_Extrema''  ] \n')
    fprintf('              N: [ possitive scalar  |  {10} ] \n')
    fprintf('              T: [ possitive scalar or three component vector |  {10} ] \n')
    fprintf('          Order: [ odd integer  |  {3}  ] \n')
    fprintf('           in_N: [ possitive scalar  |  {20}  ] \n')
    fprintf('        in_Nskip: [ possitive scalar  |  {5}  ] \n')
    fprintf('       in_Order: [ odd integer  |  {3}  ] \n')
    fprintf('Differentiation:  [ ''diff_2p''  |  {''diff_5p''} ] \n')
    fprintf('           MaxN: [ possitive scalar  |  {5000}  ] \n')
    fprintf('           IMFs: [ possitive scalar  |  {50}  ] \n')
    fprintf('           Disp: [ 0  |  {1}  ] \n')
    return
elseif isa(varargin{1},'struct')
    optionsdefaults=varargin{1};
    varargin(1)=[];
    
elseif isa(varargin{1},'function_handle')
    funcname = lower(func2str(varargin{1}));
    if strcmp(funcname,'di')==1
        optionsdefaults=struct('Stopping','Nsifts_after_IMF',...
            'N',10,...
            'T',[],...
            'MaxN',5000,...
            'IMFs',50,...
            'Disp',1,...
            'Order',3,...
            'Nodes','DI_extrema',...
            'Differentiation','diff_5p',...
            'in_N',20,...
            'in_Nskip',5,...
            'in_Order',3);
    else
        optionsdefaults=optionsstandart;
    end
varargin(1)=[];
else
    optionsdefaults=optionsstandart;
end



numbofargs=length(varargin);
if numbofargs==0
    options=optionsdefaults;
elseif mod(numbofargs,2)==1
    error('Arguments must occur in name-value pairs.')
else
    options=optionsdefaults;
    for arg=1:2:numbofargs
        in_arg=varargin{arg};
        in_val=varargin{arg+1};
        if ~ischar(in_arg)
            error('Expected argument %d to be a string parameter name.', arg)
        end
        findattempt=strmatch(lower(in_arg),names);
        if isempty(findattempt)
            error('There is not a parameter called ''%s''.',in_arg)
        elseif length(findattempt) > 1
            findattempt=strmatch(lower(in_arg),names,'exact')
        end
    
        if ischar(in_val)
            in_val = (deblank(in_val));
        end
        check_arg(Names{findattempt},in_val);
        options.(Names{findattempt,:}) = in_val;
    end
end



function check_arg(arg,val)
switch arg
    case {'Stopping'}
        if ~isa(val,'char') || ~any(strcmpi(val,{'Nsifts','Nsifts_AND_IMF','Nsifts_after_IMF','Single_T','Double_T'})) 
            error('Invalid value for parameter ''%s''.',arg);
        end
        
    case {'N','MaxN','IMFs','in_N','in_Nskip','Disp'}
        if isnumeric(val)
            if round(val)~=val
               error('Invalid value. Parameter %s takes integer inputs.',arg)
            end
        else
            error('Invalid value. Parameter %s takes integer inputs.',arg)
        end

    case {'Order','in_Order'}
        if isnumeric(val)
            if mod(val,2)==0
               error('Invalid value. Parameter %s takes odd numbers.',arg)
            end
        else
            error('Invalid value. Parameter %s takes odd numbers.',arg)
        end

    case {'T'}
        if ~isnumeric(val)
            error('Invalid value. Parameter %s takes real numbers.',arg)
        end

    case {'Differentiation'}
        if ~isa(val,'char') || ~any(strcmpi(val,{'diff_2p','diff_5p'})) 
            error('Invalid value for parameter ''%s''.',arg);
        end

    otherwise
        error('Unknown parameter name.')
end    



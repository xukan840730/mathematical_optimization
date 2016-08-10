function [out coeffnames variables] = mapcoeffstoports(this,varargin)
%MAPCOEFFSTOPORTS

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/05/23 08:14:38 $

coeffnames = [];
variables = [];
out = parse_mapcoeffstoports(this,varargin{:});

if strcmpi(out,'on')
    nstages = length(this.Stage);
    
    % User-defined coefficient names
    idx = find(strcmpi(varargin,'CoeffNames'));
    if ~isempty(idx)&&~isempty(varargin{idx+1})
        coeffnames = varargin{idx+1};
        
        % check if the specified coefficient name is a structure and the number
        % of coefficient names' stages match with filter stages.
        errorflag = false;
        if ~isstruct(coeffnames);
            errorflag = true;
        else
            fd = fields(coeffnames);
            if length(fd)~=nstages, errorflag=true; end;
        end
        
        if errorflag
            error(generatemsgid('InvalidValue'), ...
                ['The CoeffNames value of the realized multistage filter must '...
                 'be a structure array containing %d fields.'],nstages);
        end
        
        % check if the field name is in the format: Name.Stage1, Name.Stage2, etc.
        errorflag = false;
        for stage = 1:length(fd)
            stagename = fd{stage};
            idx = findstr(stagename,'Stage');   % Find 'Stage' string
            if ~isempty(idx)
                % check that the stage number is appended after the string
                % 'Stage' i.e. Stage1, Stage2.
                stagenum = str2double(stagename(idx+5:end));
                if ~isfinite(stagenum)||(stagenum<1)||(stagenum>nstages)
                    errorflag = true;
                end
            else
                errorflag = true;
            end
        end
        % provide error message when the name is not Stage1, Stage2, etc.
        if errorflag
            error(generatemsgid('InvalidValue'), ...
                ['The CoeffNames value of a multistage filter must be a structure with fields ',...
                 '''Stage<N>'' where N is the index of the filter stage. For example, a 2-stage filter ',...
                 'requires a structure with fieldnames ''Stage1'' and ''Stage2''.']);
        end
    end
    
    % Get default coefficient names and coefficient variables.
    for k=1:nstages
        [temp stagecoeff stagevar] = mapcoeffstoports(this.Stage(k),'MapCoeffsToPorts','on');
        coeffs.(sprintf('Stage%d',k))= stagecoeff;
        variables.(sprintf('Stage%d',k)) = stagevar;
    end
    
    % If no user-specified coefficient names, use default.
    if isempty(coeffnames)
        coeffnames = coeffs;
    end
    
end

% [EOF]

function [vals, errStr] = evaluatevars(strs, name)
%EVALUATEVARS   Evaluate variables in the MATLAB workspace.
%
%   EVALUATEVARS will take a string (or cell array of strings)
%   representing filter coefficients or the Sampling frequency and evaluate 
%   it in the base MATLAB workspace. If the variables exist and are numeric, the 
%   workspace variables values are returned in VALS, if they do not exist, an 
%   error dialog is launched and the error message is returned in ERRSTR. 
%
%   Input:
%     strs   - String or cell array of strings from edit boxes
%     names  - String or cell array of names for the edit boxes.  This
%              allows EVALUATEVARS to give customized error messages if the
%              editboxes are empty.  If this input is not given a generic
%              message 'Editboxes cannot be empty.' will be given.  If this
%              input is empty it will be ignored.
%
%   Outputs:
%     vals   - Values returned after evaluating the input strs in the
%              MATLAB workspace.
%     errStr - Error string returned if evaluation failed.

%   Author(s): R. Losada, P. Costa
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.9.4.6 $  $Date: 2009/07/27 20:31:55 $ 

errStr = '';
vals = {};

if  iscell(strs)
    for n = 1:length(strs), % Loop through strings
        if ~isempty(strs{n}),
            try
                vals{n} = evalin('base',['[',strs{n},']']);
                % Check that vals is a numeric array and not a string.
                if ~isnumeric(vals{n}), errStr=[strs{n} ' is not numeric.']; end     
            catch
                errStr = sprintf('The variable %s is not defined in the MATLAB workspace.', strs{n});
                break;
            end
        else
            if nargin > 1 & ~isempty(name{n})
                errStr = sprintf('The %s edit box cannot be empty.', name{n});
            else
                errStr = ['Edit boxes cannot be empty.'];
            end
            break;
        end 
    end
else
    if ~isempty(strs),
        try
            vals = evalin('base',['[',strs,']']);
            if ~isnumeric(vals), errStr=[strs ' is not numeric.']; end
        catch
            errStr = sprintf('The variable %s is not defined in the MATLAB workspace.', strs);
        end
    else
        if nargin > 1 & ~isempty(name)
            errStr = sprintf('The %s edit box cannot be empty.', name);
        else
            errStr = ['Edit boxes cannot be empty.'];
        end
    end
end

if nargout < 2,
    if ~isempty(errStr), error(generatemsgid('SigErr'),errStr); end % Top level try catch will display the error dialog.
end

% [EOF] evaluatevars.m

function variables = checkoptimizescalevalues(this,variables)
%CHECKOPTIMIZESCALEVALUES check if optimize scale values is possible

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/07/14 04:00:37 $


issvnoteq2one = this.issvnoteq2one;
if this.OptimizeScaleValues && ~all(issvnoteq2one),
    % Unit scale values cannot be skipped when specified through a port
    warning(generatemsgid('UnitScaleValues'), ...
        'Unable to optimize unit scale values when specified through a port.');
    that = copy(this);
    that.OptimizeScaleValues = false;
    g = that.privScaleValues.';
    variables{3} = g;
end


% [EOF]

function fr = whichframes(this)
%WHICHFRAMES   Returns the frames for the classic iir.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2.54.1 $  $Date: 2009/12/17 13:59:25 $

fr = dmom_whichframes(this);

if ~isspecify(this)
    fr(end).constructor = 'siggui.ellipoptsframe';
    fr(end).setops      = {'MatchExactly', get(this, 'MatchExactly')};
end

% [EOF]

function str = gettitlestr(this)
%GETTITLESTR   PreGet function for the 'titlestr' property.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2.54.1 $  $Date: 2009/12/17 13:59:20 $

str = sprintf('%% %s %s filter designed using FDESIGN.%s.', ...
    get(this, 'Tag'), ...
    get(this, 'ResponseType'), ...
    upper(get(this, 'ResponseType')));

% [EOF]

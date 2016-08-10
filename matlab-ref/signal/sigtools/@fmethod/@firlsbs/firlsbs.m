function this = firlsbs
%FIRLSBS   Construct a FIRLSBS object.

%   Author(s): J. Schickler
%   Copyright 1999-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2006/06/27 23:40:04 $

this = fmethod.firlsbs;

set(this, 'DesignAlgorithm', 'FIR least-squares');

% [EOF]

function str = designdesc(d)
%DESIGNDESC   Returns the design comment.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.2.4.2.54.1 $  $Date: 2009/12/17 13:59:19 $

str = sprintf('%% Construct an FDESIGN object and call its %s method.', ...
    upper(designfunction(d)));

% [EOF]

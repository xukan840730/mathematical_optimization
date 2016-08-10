function g = thisnormalize(Hd)
%THISNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/04/12 23:58:24 $

lat = Hd.reflattice;
lad = Hd.refladder;
g = [max(abs(lat)) max(abs(lad))];
Hd.reflattice = lat/g(1);
Hd.refladder = lad/g(2);


% [EOF]

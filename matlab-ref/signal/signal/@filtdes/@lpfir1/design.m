function  Hd = design(h,d)
%DESIGN Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2002/04/15 00:25:34 $

% Set up design params
N = get(d,'order');
Fc = get(d,'Fc');
win = generatewindow(d);

scaleflag = determinescaleflag(d);

b = fir1(N,Fc,'low',win,scaleflag);

% Construct object
Hd = dfilt.dffir(b);



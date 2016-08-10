function [y,zf] = df2sosfilter(q,num,den,sv,issvnoteq2one,x,zi)
% DF2SOSFILTER Filter for DFILT.DF2SOS class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2004/04/12 23:57:26 $

x = quantizeinput(q,x);

% Initialize
y = x;
zf = zi;

% Form 1/a0
den(:,1) = den(:,1).^-1;

% Vectorize num and den (section by section)
num = num.';
num = num(:);
den = den.';
den = den(:);

% Make sure that sv is not empty.  It doesn't matter what we set it to,
% since we only multiply the sv value based on the issvnoteq2one vector.
% This was causing problems when allocating a STL complex object for
% scalevalues (see: mat2cpp in
% toolbox/filterdesign/quanitzation/filters/include/flutilities.h)
if isempty(sv), sv = 0; end

% Make num and den row vectors
num = num.';
den = den.';

[y,zf] = df2sosfilter(num,den,sv,issvnoteq2one,x,zi);

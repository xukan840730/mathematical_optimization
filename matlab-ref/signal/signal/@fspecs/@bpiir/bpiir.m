function this = bpiir(varargin)
%BPIIR   Construct a BPIIR object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/06/16 08:28:51 $

this = fspecs.bpiir;

respstr = 'Bandpass';
fstart = 3;
fstop = 5;
nargsnoFs = 8;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});
% [EOF]

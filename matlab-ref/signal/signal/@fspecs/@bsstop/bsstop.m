function h = bsstop(varargin)
%BSSTOP   Construct a BSSTOP object.
%   H = BSSTOP(N,Fstop1,Fstop2,Astop,Fs) constructs a bandstop filter
%   specifications object with stopband-edge specifications.
%
%   N is the filter order and must be an even positive integer.
%
%   Fstop1 is the lower stopband-edge frequency and must be a positive
%   scalar between 0 and 1 if no sampling frequency is specified or between
%   0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fstop2 is the higher stopband-edge frequency and must be a positive
%   scalar larger than Fstop1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Astop is the minimum stopband attenuation and it must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/14 15:14:04 $

error(nargchk(0,5,nargin,'struct'));
h = fspecs.bsstop;
constructor(h,varargin{:});
h.ResponseType = 'Bandstop with stopband-edge specifications.';
% [EOF]

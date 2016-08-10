function Hd = kaiserwin(this, varargin)
%KAISERWIN   Design an FIR filter using the Kaiser window.
%   Hd = KAISERWIN(Hs) designs an FIR filter that meets the specifications
%   in Hs, using a Kaiser window.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/10/23 18:48:42 $

Hd = design(this, 'kaiserwin', varargin{:});
h = getfmethod(Hd);

if ishp(this),
    Hd = firlp2hp(Hd);
    % Reset the contained FMETHOD.
    Hd.setfmethod(h);
end

% [EOF]

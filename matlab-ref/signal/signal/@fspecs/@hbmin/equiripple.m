function Hd = equiripple(this, varargin)
%EQUIRIPPLE   Design an equiripple filter.
%   Hd = EQUIRIPPLE(Hs) designs an equiripple FIR filter that meets the
%   specifications in Hs.
%
%   Hd = EQUIRIPPLE(Hs,'Minphase',true) designs a minimum-phase equiripple
%   FIR filter.
%
%   EQUIRIPPLE(...) with no uses FVTool to visualize the filter that is
%   designed.

%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.1.8.8 $  $Date: 2007/10/23 18:48:35 $

Hd = design(this, 'equiripple', varargin{:});
h = getfmethod(Hd);

if ishp(this),
    Hd = firlp2hp(Hd);
    
    % Reset the contained FMETHOD.
    Hd.setfmethod(h);
end

% [EOF]

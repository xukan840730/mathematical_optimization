function this = phase(varargin)
%PHASE   Construct a PHASE object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/14 15:10:28 $

error(nargchk(0,8,nargin,'struct'));

% Create and initialize object.
this = dspdata.phaseresp;

set(this,'Name','Phase');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,'FrequencyUnits','Hz');
% From the help of TFESTIMATE and MSCOHERE we are deducing that there are
% no units for the magnitude:
%
% The magnitude squared coherence Cxy is given by
%         Cxy = (abs(Pxy).^2)./(Pxx.*Pyy)
set(this.Metadata,'DataUnits','');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

% [EOF]

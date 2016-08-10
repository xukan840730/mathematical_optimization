function h = plotfcn(this, fcn, varargin)
%PLOTFCN   Plot engine for all of the signal's methods.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/12/27 21:29:55 $

% Search through VARARGIN for 'Parent'.
hax = [];
if nargin > 2
    indx = find(strcmpi('parent', varargin));
    if ~isempty(indx)
        indx = max(indx);
        hax = varargin{indx+1};
        varargin(indx:indx+1) = [];
    end
end

% If we did not find an axes to use.
if isempty(hax)
    hax = newplot;
end

d  = get(this, 'Data');

if this.NormalizedFrequency
    t = [0:length(d)-1];
    xlbl = 'samples';
else
    fs = get(this, 'Fs');

    t = [0:1/fs:length(d)/fs-1/fs];
    [t, m, xlbl] = engunits(t, 'latex', 'time');
end

h = feval(fcn, t, d, 'Parent', hax, varargin{:});

ylabel(hax, xlate('Amplitude'));
xlabel(hax, xlate(sprintf('Time (%s)', xlbl)));

% [EOF]

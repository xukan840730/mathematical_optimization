function [P, W] = noisepsd(this, L, opts, optsstruct)
%NOISEPSD   Calculate the noise psd.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2009/07/27 20:29:10 $

if nargin < 4
    optsstruct.showref = true;
    optsstruct.sosview = [];
    if nargin < 3
        opts = dspopts.spectrum;
        opts.NFFT = 512; % NFFT default is 'Nextpow2', but here we want a numeric value.
        if nargin < 2
            L = 10;
        end
    end
end

% If there is more than 1 filter, we ignore the sosview settings.
if length(this) > 1 || ~isa(this(1).Filter, 'dfilt.abstractsos')
    optsstruct.sosview = [];
end

fs = getmaxfs(this);
if isempty(fs)
    fs = 2*pi;
end

if ishalfnyqinterval(opts);
    wmin = 0;
    wmax = fs/2;
elseif opts.CenterDC
    wmin = -fs/2;
    wmax = fs/2;
else
    wmin = 0;
    wmax = fs;
end

% Always get the data from normalized and twosided.  COMPLETEFREQRESP will
% take care of the rest.
opts.SpectrumType        = 'twosided';
opts.NormalizedFrequency = true;
opts.CenterDC            = false;

nfft = opts.NFFT;

for indx = 1:length(this)
    hindx = this(indx).Filter;
    
    fs = get(this(indx), 'Fs');
    if isempty(fs), fs = getmaxfs(this); end
    if isempty(fs), fs = 2*pi;           end

    if ~isempty(optsstruct.sosview)
        hindx = getfilters(optsstruct.sosview, hindx);
    end
    
    if any(isquantized(hindx)) && optsstruct.showref
        hindx = [hindx reffilter(hindx)];
    end
    for jndx = 1:length(hindx)
        opts.NFFT = force2even(max(4, round(nfft*fs/(wmax-wmin))));
        Hpsd = noisepsd(hindx(jndx), L, opts);
        P{indx}(:, jndx) = Hpsd.Data;
    end
    P{indx} = convert2db(P{indx})/2;

    % Complete the response for the filter from wmin to wmax based on
    % the sampling frequency of this individual filter.
    [P{indx}, W{indx}] = completefreqresp(P{indx}, fs, wmin, wmax);
end

% --------------------------------------------------------
function nfft = force2even(nfft)

if rem(nfft, 2)
    nfft = nfft+1;
end

% [EOF]

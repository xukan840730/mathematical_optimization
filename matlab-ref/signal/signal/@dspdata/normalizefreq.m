function normalizefreq(this)
%NORMALIZEFREQ   Normalize frequency specifications.
%   NORMALIZEFREQ(H) will normalize the frequency specifications of the DSP
%   data object (<a href="matlab:help dspdata">dspdata</a>) H. This sets the data object's
%   'NormalizedFrequency' property to true, and normalizes the Frequencies
%   by Fs, so that they lie between 0 and 1.
%
%   NORMALIZEFREQ(H,BFL) sets the 'NormalizedFrequency' property to the
%   value of the boolean flag, BFL, which can be either true or false. If
%   not specified, BFL defaults to true. If BFL is false, the frequencies
%   in the data object H will be converted to linear frequencies.
%
%   NORMALIZEFREQ(H,false,Fs) allows for the setting of a new sampling
%   frequency, Fs, when the 'NormalizedFrequency' property is set to false.
%
%   EXAMPLE: Use the periodogram to estimate the power spectral density of
%            % a noisy sinusoidal signal. Then store the results in a PSD
%            % data object, convert it to normalized frequency, and plot it.
%  
%            Fs = 32e3;   t = 0:1/Fs:2.96;
%            x = cos(2*pi*t*1.24e3)+ cos(2*pi*t*10e3)+ randn(size(t));
%            Pxx = periodogram(x);
%            hpsd = dspdata.psd(Pxx,'Fs',Fs) % Create a PSD data object.
%            normalizefreq(hpsd);            % Normalize frequencies.
%            plot(hpsd)
%
%   See also DSPDATA, SPECTRUM.

%   Author(s): P. Pacheco
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2008/09/13 07:14:45 $

% Help for NORMAZLIZEFREQ method.

% [EOF]

function twosided(this)
%TWOSIDED   Convert a one-sided spectrum to a two-sided spectrum.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/10/18 21:05:25 $

newSpectrumType = 'twosided';
if strcmpi(this.SpectrumType,newSpectrumType),
    return;    % Spectrum already two-sided.
end

if this.NormalizedFrequency,
    fnyq = pi; 
else
    fnyq = this.getfs/2;
end

Pxx = this.Data;
W   = this.Frequencies;
[Nfft,nchans] = size(Pxx);

% Rebuild the 'twosided' PSD from the 'onesided' PSD.
startIdx = Nfft+1;
if isevenwholenfft(this,Nfft,W),      % EVEN "whole" NFFT
    endIdx = (Nfft-1)*2;
    Pxx(2:Nfft-1,:) = Pxx(2:Nfft-1,:)/2;
    Pxx(startIdx:endIdx,:) = Pxx(Nfft-1:-1:2,:);  % Add positive half.
    W(startIdx:endIdx) = W(2:Nfft-1)+fnyq;
    
else                                  % ODD "whole" NFFT
    endIdx = (Nfft*2)-1;
    Pxx(2:Nfft,:) = Pxx(2:Nfft,:)/2;
    Pxx(startIdx:endIdx,:) = Pxx(Nfft:-1:2,:);    % Add positive half.
    W(startIdx:endIdx) = W(2:Nfft)+fnyq;
end

this.Data = Pxx;
this.Frequencies = W;
setspectrumtype(this,newSpectrumType); % Uses priv property to produce better error msg.

% [EOF]

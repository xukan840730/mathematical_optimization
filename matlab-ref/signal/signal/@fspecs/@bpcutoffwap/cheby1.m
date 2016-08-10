function Hd = cheby1(this, varargin)
%CHEBY1 Chebyshev Type I digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/16 08:28:20 $

N  = this.FilterOrder;
Ap = this.Apass;

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

Fc1 = this.F3dB1;
Fc2 = this.F3dB2;

normalizefreq(this, nfreq);

if Fc2<=Fc1,
    error(generatemsgid('invalidSpec'), 'The frequency specifications {F3dB1, F3dB2} must have increasing values.');
end

% Make up a 3dB point for digital lowpass (arbitrarily)
Fc = (Fc1+Fc2)/2; 

% Compute theoretical lowpass passband edge
Wc = tan(pi*Fc/2);

% Find corresponding analog passband-edge frequency
Wp = Wc/cosh(1/(N/2)*acosh(1/sqrt(10^(Ap/10)-1))); % Use half the order

% Convert analog passband-edge frequency to digital
Fp = 2*atan(Wp)/pi;

% Compute alpha
alpha = cos((Fc2+Fc1)*pi/2)/cos((Fc2-Fc1)*pi/2); 
% Compute k
k = tan(Fc*pi/2)/tan((Fc2-Fc1)*pi/2);

c1 = 2*alpha*k/(k+1);
c2 = (k-1)/(k+1);

% Solve LP to BP inverse Constantinides (Antoniou p.243)
j = complex(0,1);
z = exp(j*Fp*pi);
b = c1*(z+1);
ax2 = 0.5/(c2*z+1);
b2_4ac = sqrt(c1^2*(1+z)^2-4*(c2*z+1)*(z+c2));
Fp1 = abs(real(log(ax2*(b + b2_4ac))/(j*pi)));
Fp2 = abs(real(log(ax2*(b - b2_4ac))/(j*pi)));

% Convert to bandpass with passband-edge specifications
hs = fspecs.bppass(N,Fp1,Fp2,Ap);

Hd = cheby1(hs,varargin{:});

% [EOF]

function [sos,g,Ws] = alpastop(h,N,Wp,Apass,Astop)
%ALPASTOP   

%   Author(s): R. Losada
%   Copyright 1999-2005 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2005/06/16 08:42:28 $

% Compute q, k from Apass, Astop

D = (10^(0.1*Astop)-1)/(10^(0.1*Apass)-1);
[q,k] = computeq2(h,N,D);

% Find Ws
Ws = Wp/k; % Can return as a measurement

Wc=sqrt(Wp*Ws);

% Design prototype
[sos,g] = apspecord(h,N,Wp/Wc,Apass,k,q); 

% Make transformation s -> s/Wc
sos = stosbywc(h,sos,Wc);

% [EOF]

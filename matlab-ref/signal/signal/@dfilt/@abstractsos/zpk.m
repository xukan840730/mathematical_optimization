function [z,p,k] = zpk(Hd)
%ZPK  Discrete-time filter zero-pole-gain conversion.
%   [Z,P,K] = ZP(Hd) returns the zeros, poles, and gain corresponding to the
%   discrete-time filter Hd in vectors Z, P, and scalar K respectively.
%
%   See also DFILT.   
  
%   Author: R. Losada, J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
%   $Revision: 1.2.4.2 $  $Date: 2004/04/12 23:52:45 $

[msgid,msg] = warnsv(Hd);
if ~isempty(msg),
    warning(msgid,msg);
end

sosm = get(Hd, 'sosMatrix');

z = [];
p = [];
k = prod(Hd.ScaleValues);
for indx = 1:size(sosm,1)
  [z1,p1,k1] = sos2zp(sosm(indx,:));
  z = [z;z1];
  p = [p;p1];
  k = k*k1;
end

% [EOF]

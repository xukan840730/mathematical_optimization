function [NL, NextIPorts, NextOPorts, mainparams]=delayheadconnect(q,NL,H,mainparams)
%DELAYHEADCONNECT specifies connection and quantization parameters in the
%conceptual head stage

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/05/23 08:14:27 $

% specify the connection
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);

% specify interstage connection
NextIPorts=[];
NextOPorts=[];

% [EOF]

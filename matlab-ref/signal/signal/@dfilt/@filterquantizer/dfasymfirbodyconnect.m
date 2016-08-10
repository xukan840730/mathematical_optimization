function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=dfasymfirbodyconnect(q,NL,H,mainparams)
%DFASYMFIRBODYCONNECT specifies the connection and quantization parameters in the
%conceptual body stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/12/26 22:06:16 $

%sum
set(NL.nodes(2),'qparam','double');
set(NL.nodes(5),'qparam','double');

%gain
set(NL.nodes(3),'qparam','double');

%make the connection
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(4,1,2,2);
NL.connect(3,1,5,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(5,2)];
PrevOPorts = [filtgraph.nodeport(4,1)];
NextIPorts = [filtgraph.nodeport(4,1)];
NextOPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(5,1)];

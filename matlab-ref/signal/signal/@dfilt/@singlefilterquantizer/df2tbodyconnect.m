function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=df2tbodyconnect(q,NL,H,mainparams)
%DF2TBODYCONNECT specifies the connection and quantization parameters in the
%conceptual body stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/07/27 20:29:59 $

%sum
set(NL.nodes(2),'qparam','single');
set(NL.nodes(3),'qparam','single');

%gain
set(NL.nodes(1),'qparam','single');
set(NL.nodes(4),'qparam','single');

%make the connection
NL.connect(6,1,1,1);
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(4,1,3,2);
NL.connect(3,1,5,1);
NL.connect(7,1,4,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(6,1) filtgraph.nodeport(7,1)];
PrevOPorts = [filtgraph.nodeport(5,1)];
NextIPorts = [filtgraph.nodeport(2,2)];
NextOPorts = [filtgraph.nodeport(6,1) filtgraph.nodeport(7,1)];

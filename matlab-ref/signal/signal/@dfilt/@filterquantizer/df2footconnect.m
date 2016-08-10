function [NL, PrevIPorts, PrevOPorts, mainparams]=df2footconnect(q,NL,H,mainparams);
%DF2FOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/12/26 22:05:58 $

%sum
set(NL.nodes(1),'qparam','double');
set(NL.nodes(5),'qparam','double');

%gain
set(NL.nodes(2),'qparam','double');
set(NL.nodes(4),'qparam','double');

NL.connect(2,1,1,2);
NL.connect(3,1,2,1);
NL.connect(3,1,4,1);
NL.connect(4,1,5,2);
NL.connect(5,1,6,1);

% specify the interstage connection
% last layer, therefore no next input and output 
PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,1) filtgraph.nodeport(5,1)];
PrevOPorts = [filtgraph.nodeport(1,1)];

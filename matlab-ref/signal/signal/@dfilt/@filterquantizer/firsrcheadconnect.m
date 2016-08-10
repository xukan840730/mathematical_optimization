function [NL, NextIPorts, NextOPorts, mainparams] = firsrcheadconnect(q,NL,H,mainparams,interp_order,flag)
%FIRSRCHEADCONNECT 

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2008/05/31 23:26:25 $


for m=1:interp_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','double');
end

inputidx = interp_order+2;
zohidx = interp_order+1;

% connect input to zoh
NL.connect(inputidx,1,zohidx,1);

% connect zoh to gains
for m=1:interp_order
    gainidx = m;
    NL.connect(zohidx,1,gainidx,1);
end

NextIPorts=[];

% flag = 1 implies that there is at least one layer betweeen header and
% output. Hence input node outputs to next layer
if flag == 1
    NextOPorts= filtgraph.nodeport(inputidx,1);
else
    NextOPorts = [];    
end

for m=1:interp_order
    NextOPorts = [NextOPorts, filtgraph.nodeport(m,1)];
end


% [EOF]

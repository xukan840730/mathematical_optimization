function [targs, strs] = convertchoices(this)
%CONVERTCHOICES

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2009/07/27 20:31:33 $

% Common filter structures for both dfilt 
strs = {'Direct-Form I',...
    'Direct-Form II',...
    'Direct-Form I Transposed',...
    'Direct-Form II Transposed',...
    'State-Space',...
    'Lattice Autoregressive Moving-Average (ARMA)'};

targs = {'df1','df2','df1t','df2t','statespace','latticearma'};

refHd = reffilter(this);

if issos(refHd),
    for indx = 1:4
        strs{indx}  = [strs{indx} ', SOS'];
        targs{indx} = [targs{indx} 'sos'];
    end
end

% FIR case
if isfir(refHd),
    [b,a] = tf(refHd);
    if a(1)==1,
        strs = {strs{5},...
            'Direct-Form FIR',...
            'Direct-Form FIR Transposed'};
        targs = {targs{5}, 'dffir', 'dffirt'};
    else
        strs = {strs{1:5},...
            'Direct-Form FIR',...
            'Direct-Form FIR Transposed'};
        targs = {targs{1:5}, 'dffir', 'dffirt'};
    end
    if isminphase(refHd),
        strs = {strs{:},'Lattice Moving-Average Minimum Phase'};
        targs = {targs{:}, 'latticemamin'};
    elseif ismaxphase(refHd),
        strs = {strs{:},'Lattice Moving-Average Maximum Phase'};
        targs = {targs{:}, 'latticemamax'};
    end
    if islinphase(refHd) && isreal(refHd),
        ftype = firtype(refHd);
        switch ftype,
            case {1,2},
                strs = {strs{:},'Direct-Form Symmetric FIR'};
                targs = {targs{:}, 'dfsymfir'};
            case {3,4},
                strs = {strs{:},'Direct-Form Antisymmetric FIR'};
                targs = {targs{:}, 'dfasymfir'};
        end
    end

else % IIR case

    % FD Tbx required for coupled-allpass conversions
    if isfdtbxinstalled,
        strs = {strs{:},...
            'Coupled-Allpass (CA) Lattice',...
            'Coupled-Allpass (CA) Lattice with Power-Complementary (PC) Output'};
        targs = {targs{:}, 'calattice', 'calatticepc'};
    end
    if isallpass(refHd),
        strs = {strs{:},'Lattice allpass'};
        targs = {targs{:}, 'latticeallpass'};
    end
end


% [EOF]

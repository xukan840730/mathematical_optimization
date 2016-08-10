function lattice = setlattice(this, lattice)
%SETLATTICE Overloaded set on the Lattice property.
  
%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.3.4.3 $  $Date: 2004/04/12 23:58:35 $

% Always make coeffs a row
lattice = set_coeffs(this, lattice);

% Store lattice as reference and check data type
set(this,'reflattice',lattice);

ncoeffs = this.ncoeffs;
oldlength=0;
if ~isempty(ncoeffs), oldlength = ncoeffs(1); end
newlength = length(lattice);
ncoeffs(1) = newlength;
this.ncoeffs = ncoeffs;

% Update the ncoeffs property of the plugin.
set_ncoeffs(this.filterquantizer, ncoeffs);

if isempty(oldlength) || newlength~=oldlength,
    reset(this);
end

% Check to see if we are loading a bad latticema.
if oldlength ~= length(this.Lattice)
    if length(this.HiddenStates) == oldlength+1
        warning(generatemsgid('corruptMATFile'), ...
            sprintf('%s\n%s', ...
                'The MAT-file you are loading appears to contain a DFILT saved', ...
                'in a previous version of MATLAB.  Please resave the MAT-file.'));
        this.HiddenStates = [this.HiddenStates; 0];
    end
end

% Quantize the coefficients
quantizecoeffs(this);

% Hold an empty to not duplicate storage
lattice = [];

% [EOF]

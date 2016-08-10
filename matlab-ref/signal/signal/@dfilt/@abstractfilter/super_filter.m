function y = super_filter(Hd,x,dim)
%FILTER Discrete-time filter.
%   Y = FILTER(Hd,X) filters the data X using the discrete-time filter Hd.
%
%   Y = FILTER(Hd,X,DIM) filters array X along dimension DIM.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.7.4.8 $  $Date: 2007/12/14 15:07:00 $

error(nargchk(1,3,nargin,'struct'));

if nargin<2, x = []; end
if nargin<3, dim = []; end

if isempty(x), 
  y = x;
  return; 
end

s = size(x);
[x,perm,nshifts] = shiftdata(x,dim);
s_shift = size(x); % New size
x = reshape(x,size(x,1),[]); % Force into 2-D

% At this point, x is a 2-D matrix and we always filter along the columns
[Mx,Nx] = size(x);
if log2(Mx*Nx)>31, 
    error(generatemsgid('InvalidInput'), 'The input of the filter must contain at most 2^31 elements.');
end

nchannels = Hd.nchannels;

if ~Hd.PersistentMemory,
    % Reset the filter
    reset(Hd);
else
	if ~isempty(nchannels) && Nx ~= nchannels
		error(generatemsgid('InvalidDimensions'),'The number of channels cannot change when ''PersistentMemory'' is ''true''.');
	end
end

% Set number of channels
Hd.nchannels = Nx;


zi = Hd.HiddenStates;
% Expand the states for the multichannel case
zi = ziexpand(Hd,x,zi);

xi = getnonprocessedsamples(Hd);
BL = blocklength(Hd);

[Mxi,Nxi] = size(xi);
rs = mod(Mx+Mxi,BL); % Remaining samples
    
[y,zf] = secfilter(Hd,[xi;x(1:end-rs,:)],zi);
ly = size(y,1);
Hd.NumSamplesProcessed = Hd.NumSamplesProcessed+prod(Nx*(Mx+Mxi-rs));
setnonprocessedsamples(Hd,x(end-rs+1:end,:));

Hd.HiddenStates = zf;

if isempty(dim),
    dim = find(s>1,1,'first');
    if isempty(dim), dim = 1; end
end
s(dim) = ly;
s_shift(1) = ly;
y = reshape(y,s_shift); % Back to N-D array
y = unshiftdata(y,perm,nshifts);

y = reshape(y,s);

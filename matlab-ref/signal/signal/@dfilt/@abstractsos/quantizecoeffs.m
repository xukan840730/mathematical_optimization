function quantizecoeffs(h,eventData)
% Quantize coefficients


%   Author(s): R. Losada
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2008/05/12 21:35:43 $

if isempty(h.refsosMatrix) || isempty(h.refScaleValues)
    return;
end

% Quantize the coefficients
if strcmpi(class(h.filterquantizer), 'quantum.singlefilterquantizer'),
    if any(h.refsosMatrix(:,4)~=1),
        error(generatemsgid('invalidsosMatrix'), ...
            ['The 4th column of the sosMatrix must contain only 1.']);
    end
end
[numq,denq,svq] = quantizecoeffs(h.filterquantizer,h.refsosMatrix(:,1:3), ...
    h.refsosMatrix(:,4:6),h.refScaleValues,h.OptimizeScaleValues);
h.privNum = numq;
h.privDen = denq;
isnoteq2one = svnoteq2one(h.filterquantizer, h.refScaleValues);

if (h.OptimizeScaleValues == 1)
    if ~isempty(isnoteq2one),
        svq = svq(isnoteq2one); % Skip scale values that are equal to one
    end
else
    isnoteq2one = true(size(isnoteq2one));
end
    
h.privScaleValues = svq; 
h.issvnoteq2one = isnoteq2one;


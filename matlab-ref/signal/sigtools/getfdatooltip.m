function tip = getfdatooltip(index)
%GETFDATOOLTIP   Return a tip for FDATool given an index.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/08/03 21:40:03 $

% Define the valid for the signal only case.
validtips = [0 6 12 14 16 17 18 19 20 21 23 25];

% Add the plug-in tips.
if isfdtbxinstalled
    validtips = [validtips 1 2 4 5 8 9 10 11 13 24];
end
if issimulinkinstalled
    validtips = [validtips 15];
end
if isfdhdlcinstalled
    validtips = [validtips 3 22];
end
if isccslinkinstalled
    validtips = [validtips 7];
end

% Sort them so they show in the correct order.
validtips = sort(validtips);

ntips = length(validtips);

% Index into the VALIDTIPS vector to find the next tip.
switch validtips(rem(index+ntips-1, ntips)+1)
    case 0
        tip = sprintf('%s %s', ...
            'You can generate M-code functions of your design by going', ...
            'to the ''File'' menu and selecting ''Generate M-file''.');
    case 1
        tip = sprintf('%s %s %s %s %s', ...
            'You can quantize filters you', ...
            'have designed or imported to either fixed-point arithmetic or', ...
            'single-precision floating-point arithmetic.  Click the', ...
            '''Set Quantization Parameters'' sidebar button.  (Note:', ...
            'fixed-point arithmetic also requires the Fixed-Point Toolbox.');
    case 2
        tip = sprintf('%s %s %s', ...
            'If you have quantized a filter,', ...
            'its reference response is shown.  To disable this, go to the', ...
            '''View'' menu and uncheck ''Show Reference Filters''.');
    case 3
        tip = sprintf('%s %s  %s', 'You can generate synthesizable VHDL and Verilog', ...
            'code along with test benches from fixed-point filters.', ...
            'Go to the ''Targets'' menu and select ''Generate HDL''.');
    case 4
        tip = sprintf('%s %s %s', ...
            'You can see different views of second-order-section (SOS)', ...
            'filters.  You can view the sections individually or cumulatively.', ...
            'Go to the ''View'' menu and select ''SOS View Options''.');
    case 5
        tip = sprintf('%s %s %s', ...
            'You can scale the coefficients of second-order-section (SOS)', ...
            'filters.  After designing an SOS filter, go to the ''Edit''', ...
            'menu and select ''Reorder and Scale Second-Order-Sections''.');
    case 6
        tip = sprintf('%s %s', ...
            'You can generate C-header files containing your filter coefficients by', ...
            'going to the ''Targets'' menu and selecting ''Generate C header''.');
    case 7
        tip = sprintf('%s  %s', ...
            'You can export C-code directly to the Code Composer Studio.', ...
            'Go to the ''Targets'' menu and select ''Code Composer Studio (tm) IDE''.');
    case 8
        tip = sprintf('%s  %s', ...
            'You can design decimators, interpolators and sample-rate converters.', ...
            'Click the ''Create a multirate filter'' sidebar button.');
    case 9
        tip = sprintf('%s %s', ...
            'You can view each polyphase response of multirate filters.  To', ...
            'enable this view go to the ''View'' menu and check ''Polyphase View''.');
    case 10
        tip = sprintf('%s %s %s %s', ...
            'You can design inverse-sinc filters via the design panel.  In', ...
            'the lowpass option of the Response Type select ''Inverse Sinc', ...
            'Lowpass'' to design a lowpass inverse-sinc filter.  Highpass', ...
            'inverse-sinc filters can be designed under the highpass option.');
    case 11
        tip = sprintf('%s  %s %s', ...
            'You can design raised-cosine filters via the design panel.', ...
            'In the lowpass option of the Response Type', ...
            'select ''Raised-cosine'' to design a raised-cosine filter.');
    case 12
        tip = sprintf('%s  %s %s', ...
            'You can modify filter analyses using the Analysis Parameters.', ...
            'Select an analysis from the ''Analysis'' menu, and then go to the', ...
            '''Analysis'' menu and choose ''Analysis Parameters''.');
    case 13
        tip = sprintf('%s  %s %s %s %s', ...
            'You can import and export Xilinx coefficient (.COE) files.', ...
            'To import a filter, go to the ''File'' menu and select', ...
            '''Import Filter from Xilinx Coefficient (.COE)', ...
            'File''.  To export a filter, go to the ''Targets'' menu and', ...
            'select ''Xilinx Coefficient (.COE) File''.');
    case 14
        tip = sprintf('%s  %s %s', ...
            'You can save your filters for later use with the Filter Manager.', ...
            'Click ''Store Filter ...'' to save a filter or', ...
            '''Filter Manager ...'' to launch the manager dialog.');
    case 15
        tip = sprintf('%s  %s %s', ...
            'You can export filters into Simulink(R) models.', ...
            'Go to the ''File'' menu and select ''Export to Simulink Model''', ...
            'to launch the Filter Realization pane of FDATool.');
    case 16
        tip = sprintf('%s %s', ...
            'You can change the structure of your filtering by going to the', ...
            '''Edit'' menu and selecting ''Convert Structure''.');
    case 17
        tip = sprintf('%s %s %s %s', ...
            'You can graphically edit the poles and zeros of your filter using', ...
            'the Pole/Zero Editor.  Click the ''Pole/Zero Editor'' sidebar', ...
            'button.  Special Pole/Zero Editor options are also available', ...
            'from ''Pole/Zero Editor'' on the ''Edit'' menu.');
    case 18
        tip = sprintf('%s %s', ...
            'Right clicking on the X- and Y-labels of the analyses allows', ...
            'you to change the X- or Y-units used by the plot.');
    case 19
        tip = sprintf('%s %s', ...
            'You can add data markers to your analyses by single-clicking', ...
            'on the analysis plot.');
    case 20
        tip = sprintf('%s %s %s', ...
            'You can import filters from the workspace into any supported', ...
            'filter structure.  Go to the ''File'' menu and select', ...
            '''Import Filter from Workspace''.');
    case 21
        tip = 'You can save and load an entire FDATool sessions from the ''File'' menu.';

    case 22
        tip = sprintf('%s %s %s %s', ...
            'You can generate M-code that will automate the generation of', ...
            'synthesizable VHDL and Verilog code along with test benches for', ...
            'fixed-point filters.  Go to the ''Targets'' menu and select', ...
            '''Generate HDL'' then check ''Generate M-file''.');

    case 23
        tip = sprintf('%s %s %s', ...
            'You can display filter coefficients in Decimal, Hexadecimal, and', ...
            'Binary format.  Select ''Filter Coefficients'' from the ''Analysis''', ...
            'menu, and then go to the ''Analysis'' menu and choose ''Analysis Parameters''.');

    case 24
        tip = sprintf('%s %s', ...
            'You can analyze and test your fixed-point implementations using the', ...
            '''Magnitude Response Estimate'' and ''Round-off noise power spectrum'' analyses.');
    
    case 25
        tip = sprintf('%s %s', ...
            'You can draw Spectral Masks to your magnitude plot.  Go to the', ...
            '''View'' menu and select ''User-defined Spectral Mask''.');
end

tip = xlate(tip);

% [EOF]

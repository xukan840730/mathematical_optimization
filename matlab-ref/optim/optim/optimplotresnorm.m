function stop = optimplotresnorm(x,optimValues,state)
% OPTIMPLOTRESNORM Plot value of the norm of residuals at each iteration.
%
%   STOP = OPTIMPLOTRESNORM(X,OPTIMVALUES,STATE) plots OPTIMVALUES.resnorm.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTRESNORM as the plot
%   function
%     options = optimset('PlotFcns',@optimplotresnorm);
%
%   Pass the options into an optimization problem to view the plot
%     lsqnonlin(@(x) sin(3*x),[1 4],[],[],options);

%   Copyright 2006-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2009/08/29 08:28:47 $

persistent plotavailable
stop = false;

switch state
     case 'init'
         if isfield(optimValues,'resnorm')
             plotavailable = true;
         else
             plotavailable = false;
             title(sprintf('Norm of Residuals: not available'),'interp','none');
         end
    case 'iter'
        if plotavailable
            if optimValues.iteration == 0
                % The 'iter' case is  called during the zeroth iteration,
                % but it now has values that were empty during the 'init' case
                plotresnorm = plot(optimValues.iteration,optimValues.resnorm,'kd', ...
                    'MarkerFaceColor',[1 0 1]);
                xlabel(sprintf('Iteration'),'interp','none');
                set(plotresnorm,'Tag','optimplotresnorm');
                ylabel(sprintf('Norm of residuals'),'interp','none');
                title(sprintf('Norm of Residuals: %g',norm(optimValues.resnorm)),'interp','none');
            else
                plotresnorm = findobj(get(gca,'Children'),'Tag','optimplotresnorm');
                newX = [get(plotresnorm,'Xdata') optimValues.iteration];
                newY = [get(plotresnorm,'Ydata') optimValues.resnorm];
                set(plotresnorm,'Xdata',newX, 'Ydata',newY);
                set(get(gca,'Title'),'String',sprintf('Norm of Residuals: %g',optimValues.resnorm));
            end
        end
end

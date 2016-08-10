addpath(genpath('~/matlab'))
if strcmp(strtrim(getComputerName()),'1282757')
  cd ~/Projects
else
  cd ~
end
set(0,'DefaultFigurePaperPositionMode','auto')
scnsize=get(0,'Monitorpositions');
set(0,'DefaultFigurePosition', scnsize(1,:))
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontWeight','bold');
set(0,'DefaultLineLineWidth',1.5)
colordef white


'Hello.'

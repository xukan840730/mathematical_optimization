function[bsat,isomap,pdd] = iso_calc(NCEP_bmmsat,UVic_mmsat,UVic_bmmsat)

%This routine determines whether significant periods of melting conditions
%occur anywhere in the model domain.  It does this by applying a UVic SAT 
%anomaly to a base NCEP state to obtain monthly mean SAT fields.  These
%are then analysed to see if either bounding monthly value is above 0C.  If
%neither value is above 0C, a negative isomap flag is obtained.  If one value is
%obtained, linear interpolation is carried out to determine what fraction
%of that month is above 0C.  This fraction is added to the previous
%fraction, or saved for the next monthly accumulation, to determine if the threshold
%time for melting has passed.  If so, a positive isomap flag is obtained.  If both values
%are above 0C, the entire month is used with the previous/next month's
%fraction (i.e. that month's fraction of time spent above freezing is 1.0).

%Note: attotated version of the core routine (within j loops, excluding the
%calculation of t) is in file 'isostabilis_core.m').

%initialize some values

%days per month
days_per_mon = 365.25/12.;

%HARD-CODED array sizes: assume 500 years
bsat    = zeros(13,100,100);
isomap   = zeros(500,100,100);  
pdd      = zeros(500,100,100); 
isocount = zeros(100,100);

isothresh = 1.;

%determine base state to use for interpolation.
%If at last month, bsat(2) = first month of the cycle (cyclic)
for m = 1:13;
  if (m <= 12)
    bsat(m,:,:) = NCEP_bmmsat(m,:,:) - UVic_bmmsat(m,:,:);
  elseif (m == 13)
    bsat(m,:,:) = NCEP_bmmsat(1,:,:) - UVic_bmmsat(1,:,:);
  end
end

for y = 1:499
  %zero isocount at start of each new year
  isocount(:,:) = 0.0;
  for m = 1:12;    
    for i = 1:100;
      for j = 1:50;
        %NH 'year' defined as Jan. 15-Jan. 15 to avoid truncation of melt
        %periods during NH summer.       
        t=(y-1)*12 + m; 
        sat(1) = UVic_mmsat(t,j,i)+bsat(m,j,i);
        sat(2) = UVic_mmsat(t+1,j,i)+bsat(m+1,j,i);
        
        if (sat(1) && sat(2) > 0.)
          isocount=isocount+1.;
          if (isocount >= isothresh);
            isomap(y,j,i) = 1;
          end
          pdd(y,j,i) = pdd(y,j,i) + days_per_mon*(min(sat) + 0.5*abs(sat(2)-sat(1)));
      
        elseif (sat(1) && sat(2) <= 0.)
          isocount=0.;
    
        else
          if (sat(2) > sat(1))
            isocount = sat(2)/(sat(2)-sat(1));
            if (isocount >= isothresh);
              isomap(y,j,i) = 1;
            end
            pdd(y,i,j) = pdd(y,j,i) + days_per_mon*0.5*isocount*sat(2);
          elseif (sat(1) > sat(2))
            isocount = isocount - sat(1)/(sat(2)-sat(1));
            if (isocount >= isothresh);
              isomap(y,j,i) = 1;
            end
            pdd(y,j,i) = pdd(y,j,i) + days_per_mon*0.5*isocount*sat(1);
          end
        end        
      end
      for j = 51:100;

        %SH 'year' defined as Jul. 15-Jul. 15 to avoid truncation of melt
        %periods during SH summer.
        t = (y-1)*12 + 6 + m;
        sat(1) = UVic_mmsat(t,j,i)+bsat(m,j,i);
        sat(2) = UVic_mmsat(t+1,j,i)+bsat(m+1,j,i);
        
        if (sat(1) && sat(2) > 0.)
          isocount=isocount+1.;
          if (isocount >= isothresh);
            isomap(y,j,i) = 1;
          end
          pdd(y,j,i) = pdd(y,j,i) + days_per_mon*(min(sat) + 0.5*abs(sat(2)-sat(1)));
      
        elseif (sat(1) && sat(2) <= 0.)
          isocount=0.;
    
        else
          if (sat(2) > sat(1))
            isocount = sat(2)/(sat(2)-sat(1));
            if (isocount >= isothresh);
              isomap(y,j,i) = 1;
            end
            pdd(y,i,j) = pdd(y,j,i) + days_per_mon*0.5*isocount*sat(2);
          elseif (sat(1) > sat(2))
            isocount = isocount - sat(1)/(sat(2)-sat(1));
            if (isocount >= isothresh);
              isomap(y,j,i) = 1;
            end
            pdd(y,j,i) = pdd(y,j,i) + days_per_mon*0.5*isocount*sat(1);
          end
        end        
      end
    end
  end
  y
end
           
             
           
          


    
        
          



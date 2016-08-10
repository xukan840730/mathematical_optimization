        %calculate climate at middle of month and middle of next month.
        %middle-of-month used because monthly averages are centered around
        %middle of each month.         

        sat(1) = UVic_mmsat(t,j,i)+bsat(m,j,i);
        sat(2) = UVic_mmsat(t+1,j,i)+bsat(m+1,j,i);
        
        if (sat(1) && sat(2) > 0.)
        %if entire month above zero increment isocount by whole month and
        %calculate the pdds by integrating the area above freezing for the month.
          isocount=isocount+1.;
          if (isocount >= isothresh);
            isomap(y,j,i) = 1;
          end
          pdd(y,j,i) = pdd(y,j,i) + days_per_mon*(min(sat) + 0.5*abs(sat(2)-sat(1)));
      
        elseif (sat(1) && sat(2) <= 0.)
        %if entire month below zero; reset isocount to zero, don't accumulate
        %any pdds.
          isocount=0.;
    
        else
        %determine fraction isocount of month that is above freezing, and
        %on which side of the month the period above freezing lies.  If it lies
        %on the second side of the month, begin a fresh accumulation of
        %isocount.  If it occurs on the first side of the month, add it to 
        %previous isocount from end of previous month (this must be present, 
        %based on linearity and continuity).
        %Calculate pdd by integrating the area above freezing for the month.
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
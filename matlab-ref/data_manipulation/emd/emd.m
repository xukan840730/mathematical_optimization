function [IMF] = emd(signal,options);

% Computes the IMFs using variants of EMD.
% Signal: Input signal to get decomposed.
% options: Structure having the proper parameters generated with the aid of
%         emdoptimset function. Use 'help emdoptimset', or just
%         'emdoptimset' in Matlab command line to see details.
% IMF: A Matrix containing the resulting IMFs in its lines.
% Examples:
%1) x=randn(1,1000);IMF=emd(x);
%2) options=emdoptimset('stop','Double_T','T',[0.05,0.5,0.05]);
%   x=randn(1,1000);IMF=emd(x,options);
%
% Ver. 1.0
% Yannis Kopsinis, 23/04/2008
% kopsini@ieee.org
%
% Some code (quite many actually) of the functions
% "edge_extrapolation" and "extremes" has been 'borrowed'
% from the first version of Rilling emd.m code. 
% His new emd versions with extra functionings can be found in 
% http://perso.ens-lyon.fr/patrick.flandrin/emd.html

t=1:length(signal);
if size(signal,2)==1, signal=signal.'; end
if nargin==1, options=emdoptimset; end


if strcmpi(options.Nodes,'DI_Extrema')
    if isempty(options.in_N) || isempty(options.in_Nskip) || isempty(options.in_Order) 
        error('The values of ''in_N'' and/or ''in_Nskip'' and/or ''in_Order'' are wrong.')
    end
end

if strcmpi(options.Stopping,'Single_T')
    if length(options.T)~=1
        error('The value of ''T'' is wrong.')
    end
end

if strcmpi(options.Stopping,'Double_T')
    if length(options.T)~=3
        error('The value of ''T'' is wrong.')
    end
end


NofExtr = inf; %nofextr: Number of Extrema
IMFnum = 1; %Current IMF number.


totalstopping=0;
while ~totalstopping & IMFnum<=options.IMFs
    tic
    
    iter = 0;
    h = signal; %Initial IMF value.
    Stoppingflag = 0; %Stoppinflag takes value 1 when the IMF criterion of options.stopping have been fulfilled
    zeromeaned = 0;

    while ~Stoppingflag & iter <= options.MaxN
        iter = iter + 1;

        options.iter = iter;
        siftdetails = sifting(h,t,options);
        
        if siftdetails.valid==1
            m = siftdetails.m;
            h = h - m;
        end
        
        [Stoppingflag,totalstopping,zeromeaned]=stoppingcriteria(siftdetails,h,options,zeromeaned,iter);

        % Display progress
        if options.Disp==1
            a=toc;
            if a>3
                fprintf('   sifting number: %d \n', iter)
                tic
            end
        end
    end

    if iter==1
        disp('')
    end
    
    if options.Disp==1
       fprintf('IMF number %d has been computed with %d sifting iterations\n',IMFnum,iter-1) 
    end
    
    IMF(IMFnum,:) = h;

    signal = signal - h;

    NofExtr = siftdetails.NofExtr;
    IMFnum = IMFnum + 1;
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function siftdetails=sifting(signal,t,options);

[maxi, mini, zerocross] = extremes(signal,t,options);

NofExtr = length(maxi)+length(mini);
siftdetails.NofExtr=NofExtr;
siftdetails.valid=1;
if length(maxi) >= 2 & length(mini) >= 2
    N_extrap=floor(options.Order/2)+1;
    [maxp_ex, maxv_ex, minp_ex, minv_ex] = edge_extrapolation(maxi, mini, signal, t,N_extrap);

    %The total number of extrapolated maxima and minima should be larger than the spline
    %order (It takes effect for spline orders higher than 3.)
    if length(maxp_ex) <= options.Order | length(minp_ex) <= options.Order
        m=zeros(1,length(signal));
        siftdetails.valid=0;
        return
    end
else
    m=zeros(1,length(signal));
    siftdetails.valid=0;
    return
end


envmax = interpolation(maxp_ex,maxv_ex,t,options); % envelop of the maxima
envmin = interpolation(minp_ex,minv_ex,t,options); % envelop of the minima



siftdetails.m = 1/2*(envmax + envmin);
siftdetails.envmax=envmax;
siftdetails.envmin=envmin;
siftdetails.maxi=maxi;
siftdetails.mini=mini;
siftdetails.zerocross=zerocross;
siftdetails.NofExtr=NofExtr;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxima , minima, zerocross] = extremes(signal,t,options)

if nargin < 3
    kind = 'extrema';
    differentiation='diff_2p';
else
    kind = lower(options.Nodes);
    differentiation=options.Differentiation;
end

if strcmpi(lower(kind),'di_extrema')

    N_extrap=floor(options.in_Order/2)+1;

    if mod(options.iter,options.in_Nskip)~=0 & options.iter~=1
        kind='extrema';
    else


        h=t(2)-t(1);

        dsignal=diff_5p(signal,h); %first derivative of signal.
        options.Nodes='extrema'; %It is used in the "extremes" function for the internal sifting iterations
        options.Order=options.in_Order; %It is used in the "interpolation" function for the internal sifting iterations
        for nofsift=1:options.in_N
            
            [maxifd , minifd, zcrsfd] = extremes(dsignal,t,options);
            
            if length(maxifd)>=ceil(options.in_Order/2) & length(minifd)>=ceil(options.in_Order/2)
                [maxp_ex, maxv_ex, minp_ex, minv_ex] = edge_extrapolation(maxifd, minifd, dsignal, t, N_extrap);
            else
                break % for the rare case that some extrema in the dfastsignal get lost and lead to less than the required extrema for the specific interpolation order.
            end
            envmax = interpolation(maxp_ex,maxv_ex,t,options); % envelop of the max
            envmin = interpolation(minp_ex,minv_ex,t,options); % envelop of the min

            m = 1/2*(envmax + envmin);
            dsignal=dsignal-m;
        end
        if nofsift==1
            kind='extrema';
        else
            
            dif = sign(dsignal);

            % It considers the appearence of single zeros only.
            if any(dif==0)
                zdif=find(dif==0);      %in case of zero the sign of the previous
                if zdif(1)==1, zdif(1)=[]; end
                if zdif(end)==length(dif), zdif(end)=[]; end
                dif(zdif)=dif(zdif-1);  %point is adopted.
            end

            N = length(dif);

            sdif = dif(1:N-1).*dif(2:N);
            extremespos = find(sdif==-1);
            extremestype = dif(extremespos);

            maxima = extremespos(find(extremestype>0))+1;
            minima = extremespos(find(extremestype<0))+1;
            zerocross = sort([maxifd, minifd]);
            if signal(maxima(1))<signal(minima(1))
                tmp=maxima(1);maxima(1)=minima(1);minima(1)=tmp;
            end
            if signal(maxima(end))<signal(minima(end))
                tmp=maxima(end);maxima(end)=minima(end);minima(end)=tmp;
            end
        end
    end
end




if strcmpi(kind,'extrema')
    if strcmpi(differentiation,'diff_2p')==1
         dif = diff(signal); %approximate derivative
    else
        h=t(2)-t(1);
        dif=diff_5p(signal,h);
    end        
  
    dif = sign(dif);

    % It considers the appearence of single zeros only.
    if any(dif==0)
        zdif=find(dif==0);      %in case of zero the sign of the previous
        if zdif(1)==1, zdif(1)=[]; end
        if zdif(end)==length(dif), zdif(end)=[]; end
        dif(zdif)=dif(zdif-1);  %point is adopted.
    end

    N = length(dif);

    sdif = dif(1:N-1).*dif(2:N);
    extremespos = find(sdif==-1);
    extremestype = dif(extremespos);

    maxima = extremespos(find(extremestype>0))+1;
    minima = extremespos(find(extremestype<0))+1;

    zerocross = find(signal(1:N-1).*signal(2:N)<0);




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function YI=interpolation(X,Y,XI,options)

h = diff(X);
if any(h == 0)
    [Kx,Lx]=find(h==0);
    X(Lx)=[];
    Y(Lx)=[];
end

if options.Order~=3
    order = options.Order + 1;
    kovw=order/2-1;
    knotseqmax=[X(1)*ones(1,order), X(3+kovw-1:end-1-kovw), X(end)*ones(1,order)];
    spout = spapi(knotseqmax,X,Y);
    YI = fnval(spout,XI);

elseif  options.Order==3
    YI = interp1(X,Y,XI,'spline');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stoppingflag,totalstopping,zeromeaned]=stoppingcriteria(siftdetails,h,options,zeromeaned,siftiter);

totalstopping = 0;
if siftdetails.valid==0
    crettotal = 0;
    totalstopping = 1;
elseif strcmpi(options.Stopping , 'Double_T')
    nofcross = length(siftdetails.zerocross);
    nofextr = length(siftdetails.maxi)+length(siftdetails.mini);
    a = mean(abs(siftdetails.envmax-siftdetails.envmin))/2 ;
    sigma = abs(siftdetails.m)./a;
    cret1 = mean(sigma > options.T(1)) > options.T(3);
    cret2 = any (sigma > options.T(2));
    cret3 = abs(nofextr-nofcross)>1;
    %cret3 = any(h(siftdetails.maxi)<0) | any(h(siftdetails.mini)>0);
    crettotal = ( cret1 | cret2 | cret3) ;
    
elseif strcmpi(options.Stopping , 'Single_T')
    h_before = h + siftdetails.m;
    SD=sum( (h_before - h).^2 / h_before.^2 );
    crettotal = SD > options.T;

elseif strcmpi(options.Stopping,'Nsifts')==1
    if options.N >= siftiter
        crettotal = 1;
    else
        crettotal = 0;
    end

elseif strcmpi(options.Stopping,'Nsifts_AND_IMF')==1
    if options.N < siftiter & (~any(h(siftdetails.maxi)<0) & ~any(h(siftdetails.mini)>0))   %abs(nofextr-nofcross)<=1
        crettotal = 0;
    else
        crettotal = 1;
    end

elseif strcmpi(options.Stopping,'Nsifts_after_IMF')==1

    if  any(h(siftdetails.maxi)<0) | any(h(siftdetails.mini)>0)
        zeromeaned = 0;
        crettotal = 1;
    else
        if options.N >= zeromeaned
            zeromeaned = zeromeaned + 1;
            crettotal = 1;
        else
            crettotal = 0;
        end
    end



end

stoppingflag=~crettotal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d=diff_5p(x,h)
% First derivative using the five points rule
% e.g. see http://www.sitmo.com/eq/455
% input
% x: equally spaced samples
% h: time distance between two samples
% output
% d: First derivative of x.
%
% Yannis Kopsinis,
% kopsinis@ieee.org

N=length(x);
d(1)=(-25.*x(1)+48.*x(2)-36.*x(3)+16.*x(4)- 3.*x(5))/(12.*h);
d(2)=(-3.*x(1)-10.*x(2)+18.*x(3)-6.*x(4)+ x(5))/(12.*h);

d(3:N-2) = (x(1:N-4)-8.*x(2:N-3)+8.*x(4:N-1)-x(5:N))/(12.*h);

d(N-1)=(-x(N-4)+6.*x(N-3)-18.*x(N-2)+10.*x(N-1)+3.*x(N))/(12.*h);
d(N)=(3.*x(N-4)-16.*x(N-3)+36.*x(N-2)-48.*x(N-1)+25.*x(N))/(12.*h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maxp_ex, maxv_ex, minp_ex, minv_ex] = edge_extrapolation(maxi, mini, signal, t, N_extr);

% This function contains many parts (possibly modified) from the EMD.m function made by
%
% G. Rilling, July 2002
%
% that computes EMD (Empirical Mode Decomposition) according to:
%
% N. E. Huang et al., "The empirical mode decomposition and the
% Hilbert spectrum for non-linear and non stationary time series analysis,"
% Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
%
% with variations reported in:
%
% G. Rilling, P. Flandrin and P. Gon?alv?s
% "On Empirical Mode Decomposition and its algorithms"
% IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
% NSIP-03, Grado (I), June 2003




% Extrapolates a number of N_extr extemes around the edges of the signal
%
% method: 'mirror' referes to the method of Rilling et al., "on empirical
%          mode decomposition and its algorithms",
%          sortly, new extrapolation methods will be added.
%
% maxi (mini): the index positions (1,2,...,N) of the maxima (minima)
% N_extr: The number of extrapolating points

% maxp_ex (minp_ex): the time positions of the maxima (minima) including
%                    the extrapolated points at the edges
% maxv_ex (minv_ex): the values of the maxima (minima) including
%                    the values of the extrapolated points at the edges

maxp = t(maxi); % maxp (minp): the time positions of the maxima (minima)
minp = t(mini);
maxv = signal(maxi); % maxv (minv): the values of the maximuma (minimuma)
minv = signal(mini);


indmax = maxi;
indmin = mini;
m = signal;
NBSYM = N_extr;
lx = length(signal);
%% copy from the Rilling m-file.
if indmax(1) < indmin(1)
    if m(1) > m(indmin(1))
        lmax = fliplr(indmax(2:min(end,NBSYM+1)));
        lmin = fliplr(indmin(1:min(end,NBSYM)));
        lsym = indmax(1);
    else
        lmax = fliplr(indmax(1:min(end,NBSYM)));
        lmin = [fliplr(indmin(1:min(end,NBSYM-1))),1];
        lsym = 1;
    end
else
    if m(1) < m(indmax(1))
        lmax = fliplr(indmax(1:min(end,NBSYM)));
        lmin = fliplr(indmin(2:min(end,NBSYM+1)));
        lsym = indmin(1);
    else
        lmax = [fliplr(indmax(1:min(end,NBSYM-1))),1];
        lmin = fliplr(indmin(1:min(end,NBSYM)));
        lsym = 1;
    end
end

if indmax(end) < indmin(end)
    if m(end) < m(indmax(end))
        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
        rmin = fliplr(indmin(max(end-NBSYM,1):end-1));
        rsym = indmin(end);
    else
        rmax = [lx,fliplr(indmax(max(end-NBSYM+2,1):end))];
        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
        rsym = lx;
    end
else
    if m(end) > m(indmin(end))
        rmax = fliplr(indmax(max(end-NBSYM,1):end-1));
        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
        rsym = indmax(end);
    else
        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
        rmin = [lx,fliplr(indmin(max(end-NBSYM+2,1):end))];
        rsym = lx;
    end
end

tlmin = 2*t(lsym)-t(lmin);
tlmax = 2*t(lsym)-t(lmax);
trmin = 2*t(rsym)-t(rmin);
trmax = 2*t(rsym)-t(rmax);

% in case symmetrized parts do not extend enough
if tlmin(1) > t(1) | tlmax(1) > t(1)
    if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,NBSYM)));
    else
        lmin = fliplr(indmin(1:min(end,NBSYM)));
    end
    if lsym == 1
        error('bug')
    end
    lsym = 1;
    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
end

if trmin(end) < t(lx) | trmax(end) < t(lx)
    if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
    else
        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
    end
    if rsym == lx
        error('bug')
    end
    rsym = lx;
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);
end

mlmax =m(lmax);
mlmin =m(lmin);
mrmax =m(rmax);
mrmin =m(rmin);


if length(mlmax)<NBSYM
    KK=NBSYM-length(mlmax);
    if length(mlmax)>1
        dml=mlmax(1)-mlmax(2);
        dtl=abs(tlmax(1)-tlmax(2));
    else
        dml=mlmax-maxv(1);
        dtl=abs(tlmax-maxp(1));
    end
    mlmax=[fliplr(mlmax(1)+dml*(1:KK)) mlmax];
    tlmax=[fliplr(tlmax(1)-dtl*(1:KK)) tlmax];
end


if length(mrmax)<NBSYM
    KK=NBSYM-length(mrmax);
    if length(mrmax)>1
        dmr=mrmax(end)-mrmax(end-1);
        dtr=abs(trmax(end)-trmax(end-1));
    else
        dmr=mrmax-maxv(end);
        dtr=abs(trmax-maxp(end));
    end

    mrmax=[mrmax mrmax(end)+dmr*(1:KK)];
    trmax=[trmax trmax(end)+dtr*(1:KK)];
end


if length(mlmin)<NBSYM
    KK=NBSYM-length(mlmin);
    if length(mlmin)>1
        dml=mlmin(1)-mlmin(2);
        dtl=abs(tlmin(1)-tlmin(2));
    else
        dml=mlmin-minv(1);
        dtl=abs(tlmin-minp(1));
    end
    mlmin=[fliplr(mlmin(1)+dml*(1:KK)) mlmin];
    tlmin=[fliplr(tlmin(1)-dtl*(1:KK)) tlmin];
end


if length(mrmin)<NBSYM
    KK=NBSYM-length(mrmin);
    if length(mrmin)>1
        dmr=mrmin(end)-mrmin(end-1);
        dtr=abs(trmin(end)-trmin(end-1));
    else
        dmr=mrmin-minv(end);
        dtr=abs(trmin-minp(end));
    end
    mrmin=[mrmin mrmin(end)+dmr*(1:KK)];
    trmin=[trmin trmin(end)+dtr*(1:KK)];
end
maxv_ex = [mlmax maxv mrmax];
maxp_ex = [tlmax maxp trmax];
minv_ex = [mlmin minv mrmin];
minp_ex = [tlmin minp trmin];



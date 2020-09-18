function soln=Mean_Field_Theory(p)
% Simulate a QCL (or any other laser) using the mean-field theory discussed
% in Burghoff's "Frequency-modulated combs as phase solitons" (ArXiv:2006.12397, 2020).
% This theory reduces the Maxwell-Bloch equations of a laser down to a
% single equation that can be integrated over the round trip of a cavity, which is then iterated.
% A GUI allows certain parameters to be adjusted on the fly if desired.

% Specifically, it uses a symmetric split-step method to find the F function in equation (5).
% It also produces a theoretical plot for the theoretical form of the soliton
% (equation (9), currently only valid for a Fabry-Perot cavity with either R1=1 or R2=1).
%
% Input: a parameter structure generated using NLFM_Params
% Output: a solution structure that stores the field evolution
% Example 1: s=Mean_Field_Theory(NLFM_Params('kpp',-1000))
%            Sets all parameters to their default but the dispersion, which is -1000 fs^2/mm
% Example 2: s=Mean_Field_Theory(NLFM_Params('gc',0,'useLinPwrApprox',1))
%            Disables gain curvature and makes the linear power
%            approximation, which leads to the phase version of the NLSE
%
% Contributors: David Burghoff, Levi Humbard

%% Load params, set constants
addpath('path_functions');
J=p.J; kpp=p.kpp; gammaK=p.gammaK; R1=p.R1; R2=p.R2; lm0=p.lm0; Lc=p.Lc; df=p.df; mu=p.mu; T1=p.T1; T2=p.T2; n=p.n; Lmod=p.Lmod; aw=p.aw; Gamma=p.Gamma; h=p.h; dtperTr=p.dtperTr; numTr=p.numTr; Crnt=p.Crnt; gc=p.gc; Ls=p.Ls; useN2N3=p.useN2N3; initphi=p.initphi; maxsave=p.maxsave; plotprogress=p.plotprogress; plottheory=p.plottheory; plotinterval=p.plotinterval;  
ec=1.60217662e-19;  % electron charge
c=2.99792458e8;
hbar=1.0545718e-34;
eps0=8.85418782e-12; % epsilon 0

%% Unit conversions
mu=mu*ec;                       % Dipole matrix element to C*m
J=J*1e4 / ec;                   % current is sheet density per time, not A/cm^2
kpp = 1*kpp*1e-30 / 1e-3;       % GVD to s^2/m
aw = aw*100;                    % waveguide loss to m^-1

%% Simulation params
Tr = 2*Lc/(c/n);                    % Round trip time
dt = Tr / dtperTr;                  % Time step size
dz = dt*c/n / Crnt;                 % Space step size
Ndz = round(Lc/dz); Nz=Ndz+1;       % # node spaces / # of nodes
Lc = dz*Ndz;                        % Cavity length
z = dz*[0:Ndz]';                    % Array containing node positions
Nt = round(numTr/h);                         % Total number of time steps
z2 = [z(1:end-1); Lc+z(1:end-1)];   % doubled cavity

%% Derived params
Psat=2*hbar^2/(mu^2*T1*T2);
g0=mu^2*Gamma*2*pi*c/lm0*T1*T2*J/(hbar*eps0*c*n*Lmod);
dw = 2*pi*df;
am = log(1/R1/R2)/(2*Lc);
p.derived = struct('Psat',Psat,'g0',g0,'dw',dw,'am',am);

if isequal(Ls,'lorentzian')
    Ls = @(w) 1./(1+i*w*T2);
elseif isequal(Ls,'parabolic');
    Ls = @(w) 1-i*w*T2+(i*w*T2).^2;
end
[fs,fss]=fftfreqs(length(z2),dz); kss=2*pi*fss;
maxLs = max(real(Ls(-c/n*kss)));

%% Plot setup
disptimer = tic;
if plotprogress
    fh=dfigure('Position',[2169 218 557 746]); ml=-Inf;
    dinput(fh,'pushbutton','','Stop',@buttonpushed); pushed=0;
    dinput(fh,'edit','J',num2str(J/(1e4 / ec)),@Jedit);
    dinput(fh,'edit','kpp',num2str(kpp/(1e-30 / 1e-3)),@kppedit);
    dinput(fh,'edit','gammaK',num2str(gammaK),@betaedit);
    dinput(fh,'pushbutton','','Reset',@reset);
    dinput(fh,'pushbutton','','Gain curvature?',@toggc);
    dinput(fh,'pushbutton','','Kick phase',@phasepert);
    dinput(fh,'edit','h',num2str(h),@hedit);
    dinput(fh,'pushbutton','','Theory?',@toggletheory);
    dinput(fh,'edit','Plot interval',num2str(plotinterval),@cplottime);
end

%% Convolution function to emulate numerical diffusion
% Set up a convolution function to prevent numerical instability (emulates
% numerical diffusion)
zw=1e-7;                            % behaves like gain curvature
[fs,fss]=fftfreqs(length(z2),dz); kss=2*pi*fss;
z2s = fftshift(z2); zi=find(z2s==0); z2s(1:zi-1)=z2s(1:zi-1)-2*Lc;
gfcn = ifftshift(1/zw/sqrt(2*pi)*exp(-1/2*(z2s/zw).^2));
Fgfcn = fft(gfcn);
nrm = real(ifft(Fgfcn.*fft(ones(size(z2))))); Fgfcn = Fgfcn./nrm;
cfcn = @(xi) (ifft(Fgfcn.*fft(xi)));
%     dinput(fh,'edit','zw',num2str(zw),@zwedit);

%% Calculate steady state power function
% Power functions
% P = Steady_State_Power();
% fn=gcf; figure(fn);
[X, P] = SSP(1e2);  % for increased accuracy increase the input of SSP
P0 = P(1);
%% Gain and K functions
geff = -aw+g0*(1-1/Psat*(P+2*flipud(P)));
Ig = cumsum(geff)*dz;
Ig(Ndz+1:end)=Ig(Ndz+1:end)+log(R2);        % mirror gain is log(R) (log R1 is implicit)
K = exp(Ig);


Kh  = zeros(4*Ndz,1); Kh(1:2:end)=K;
Kh(2:2:end)=interp1([1:2:4*Ndz+1]',[K;K(1)],[2:2:4*Ndz]); Kh=[Kh;Kh];     % two periods of K(u/2)
% function co = twoperiodconv(khi,yi) % Old Khat function (slower)
%     co = conv(khi,[yi;yi],'valid')*dz;
%     co = co(1:2*Ndz);
% end
% Khat = @(f) 1/(2*Lc)*1/2*twoperiodconv(flipud(Kh),f);
FKh = fft(flipud(Kh(1:4*Ndz)));
function ko = Khat(f)                   % Fourier version of Khat (5x faster for 1000 nodes)
    ko = 1/(4*Lc)*circshift(ifft(FKh.*fft([f;f])),[1,0])*dz;
    ko = ko(1:2*Ndz);
end
Km = 1/(2*Lc)*sum(K)*dz;

%% Linear and nonlinear split-step method algorithms
    function [N1,N2,N3,dFdz] = Nfcn(Fi)             % Nonlinear
        dFdz = ifft(ifftshift((i*kss-i*dw*n/c).*fftshift(fft(Fi))));
        N1 = -g0/2/Psat.*(c/n)*(abs(Fi).^2*Km + 2*Khat(abs(Fi).^2) -3*Km*P(1));
        N1 = N1 - i*gammaK*c/n*(abs(Fi).^2*Km + 2*Khat(abs(Fi).^2));
        
        if p.useLinPwrApprox == 0 
            N1 = N1 - g0/2/Psat*(c/n)^2*(2*T1+3  *T2)*Khat(conj(dFdz).*Fi);
            N1 = N1 - g0/2/Psat*(c/n)^2*(  T1+5/2*T2)*Khat(conj(Fi).*dFdz);
        else
            gm = g0/2/Psat*(c/n)^2*(T1+1/2*T2)/(4*Lc)*(P(end)-P(1))/P(1);
            phi = unwrap(angle(F));
            dphiend = round((phi(end)-phi(1))/(2*pi))*2*pi; phi = phi - dphiend/(2*Lc)*z2;
            N1 = N1 + i*gm*abs(F).^2.*(phi-mean(phi));
        end
        
        % Group-del
        N2 =    - g0/2/Psat*(c/n)^2*(T1+5/2*T2  )*(Khat(abs(Fi).^2)+Km*abs(Fi).^2);
        N3 =    - g0/2/Psat*(c/n)^2*(  T1+3/2*T2)*Km*F.^2;
    end
    function Dv = Dfcn(Fi)              % Linear
        DA = i/2*kpp*(c/n)^3*(-kss.^2);
        DA = DA + g0/2*c/n*(Ls(i*c/n*i*(kss-dw*n/c)).^gc-maxLs);
        Dv = ifft(ifftshift(exp(h*Tr*DA/2).*fftshift(fft(Fi))));
    end

%% Primary simulation code
[phia,fa,Pgc] = Analytic_Calculations();
if isequal(initphi,'fundamental')
    phiI = phia;
elseif isequal(initphi,'rand')
    phiI = randn(size(phia));
elseif isequal(initphi(1:6),'cosine')
    Nharm = str2num(initphi(7:end))
    phiI = 3*cos(Nharm*2*pi/(2*Lc)*z2);
elseif isequal(initphi(1:8),'harmonic')
    Nharm = str2num(initphi(9:end));
    phias = circshift(phia,[-Ndz+round(2*Ndz/Nharm/2)]);
    phiI= phias(round(mod([0:2*Ndz-1],2*Ndz/Nharm))+1);
end
F = sqrt(Pgc).*exp(i*phiI);                                % slow-intensity envelope

% Things to save
Nout = min(maxsave,Nt); iio=round(linspace(1,Nt,Nout));
[aF,aA,aphi,af,aS] = deal(zeros(length(F),Nout));
[ts,ashft]=deal(zeros(1,Nout));

ct = 0;
for ii=1:Nt 
    F = Dfcn(F);
    [N1,N2,N3,dFdz] = Nfcn(F);
    if useN2N3
        NF  = N1.*F + N2.*dFdz + N3.*conj(dFdz);
        dNFdz = ifft(ifftshift((i*kss-i*dw*n/c).*fftshift(fft(NF))));
        N2F = N1.*NF + N2.*dNFdz + N3.*conj(dNFdz);
        F = F + h*Tr*NF + 1/2*(h*Tr)^2*N2F;
    else 
        F = exp(h*Tr*N1).*F;
    end
    F = Dfcn(F);
%     F = F+(randn(size(F))+i*randn(size(F)));
    ct = ct+h*Tr; ts(ii)=ct;

    if toc(disptimer)>plotinterval & plotprogress==1
        phi = unwrap(angle(F));
        dphiend = round((phi(end)-phi(1))/(2*pi))*2*pi;
        phi = phi - dphiend/(2*Lc)*z2;
        
        [phia,fa,Pgc,Deff] = Analytic_Calculations();
        
        [~,ml]=max(phi/Deff);
        phia = circshift(phia,ml);
        fa = circshift(fa,ml);
        Pgc = circshift(Pgc,ml);
        
        [~,l6] = min(abs(z2-0.75*2*Lc));
        sfcn = @(x)circshift(x,[l6-ml]);
        
        a2=subplot(4,1,2); plot(z2/1e-3,sfcn(phi-mean(phi)));
        xyt('Position (mm)','Phase (rad)','Phase');
        if plottheory, hold all; plot(z2/1e-3,sfcn(phia-mean(phia))); hold off; end
        
        a1=subplot(4,1,1); plot(z2/1e-3,sfcn(abs(F).^2));
        xyt('Position (mm)','Intensity (V/m)^2',['Intensity, t=',num2str(ii*h),' T_r']);
        if plottheory, hold all; plot(z2/1e-3,sfcn(Pgc)); hold off; end
        ylim([0,2*P0]);
        
        a3=subplot(4,1,3); plot(z2/1e-3,sfcn(diff([NaN;phi])/dz*(-c/n)/(2*pi)/1e12));
        if plottheory, hold all; plot(z2/1e-3,sfcn(fa/1e12)); hold off; end
        xyt('Position (mm)','Frequency (THz)','Frequency');% pause
        
        s1 = abs(fftshift(fft(F))).^2;
        s2 = abs(fftshift(fft(sqrt(Pgc).*exp(i*phia)))).^2;
        s3 = real(Ls(-c/n*fss*2*pi));
        a4=subplot(4,1,4); plot(-c/n*fss/1e12,s1/max(s1)); title('Spectrum');
        xyt('Frequency (THz)','Intensity (a.u.)','Spectrum');
        if plottheory, hold all; plot(-c/n*fss/1e12,s2/max(s2)); hold off; end
        hold all; plot(-c/n*fss/1e12,s3/max(s3)); hold off;
        a4.YLim = [0,a4.YLim(2)];
        
        as=[a1,a2,a3,a4];
        for jj=1:length(as)
            as(jj).Position(1)=as(jj).Position(1)+.08;
            as(jj).Position(3)=as(jj).Position(3)-.05;
        end
        a1.Position=[0.22 0.7994 0.725 0.1577];
        a2.Position=[0.22 0.5683 0.725 0.1577];
        a3.Position=[0.22 0.3398 0.725 0.1577];
        drawnow; disptimer=tic;
        if pushed==1, break; end
    end
    if any(iio==ii)
        myi = find(iio==ii);
        phi = unwrap(angle(F));
        dphiend = round((phi(end)-phi(1))/(2*pi))*2*pi;
        phi = phi - dphiend/(2*Lc)*z2;
        [~,l6] = min(abs(z2-0.75*2*Lc));
        [~,ml]=max(phi/sign(kpp));
        aF   (:,myi)=F;
        aphi (:,myi)=phi;
        aA   (:,myi)=abs(F);
        af   (:,myi)=diff([phi(end);phi])/dz*(-c/n)/(2*pi)/1e12;
        aS   (:,myi)=fftshift(fft(F));
        ashft(myi)=l6-ml;
    end
end

% Construct solution
soln=struct();
sr = [1:myi];  % in case we ended early
soln.z=z2;                                                          % position
soln.P=P;                                                           % steady state power function
soln.t=ts(iio(sr));                                                 % times (s)
soln.F=aF(:,sr);                                                    % envelopes
soln.phi=aphi(:,sr);                                                % inst. phases
soln.A=aA(:,sr);                                                    % amplitudes
soln.f=af(:,sr);                                                    % inst. frequency
soln.shft=ashft(sr);                                                % shifts to move the jump to 75%
soln.S=aS(:,sr);                                                    % spectrum
soln.Sf=-c/n*fss;                                                   % spectrum frequencies

% Last 10% metrics 
[soln.pn10] = deal(zeros(size(fss)));                             
for jj=1:length(soln.Sf)
    phiS = unwrap(angle(soln.S(jj,:)));
    gr = [round(myi*9/10):myi];
    dphi = phiS-polyval(polyfit(soln.t(gr),phiS(gr),1),soln.t);
    soln.pn10(jj)=mean(abs(dphi(gr)).^2);                           % phase noise
end
soln.g10  = real(sqrt(1-soln.pn10));                                % g^(1) coherence (spectrally-resolved)
soln.aPS10 = mean(abs(soln.S(:,gr)).^2,2);                          % average power spectrum
soln.ag10 = sqrt(sum(soln.g10.^2.*soln.aPS10)./sum(soln.aPS10));    % average g^(1), weighted by power
                 

soln.zn = z;                                                        % normal cavity position (m)
soln.toNormal = @toNormal;                                          % function that returns forward and backward waves          
soln.params=p;                                                      % input parameters
function [Epi,Emi] = toNormal(Ei)
    Epi = [Ei(1:Ndz);NaN*Ei(Ndz+1)];
    Emi = [Ei(1)*NaN;flipud(Ei(Ndz+1:end))];
end
        

%% Plot button handling
    function buttonpushed(varargin)
        pushed=1;
    end
    function Jedit(varargin)
        J = str2num(varargin{1}.String)*1e4 / ec;
        g0=mu^2*Gamma*2*pi*c/lm0*T1*T2*J/(hbar*eps0*c*n*Lmod);
        tic; [X, P] = SSP(1e2); toc % for increased accuracy increase the input of SSP
        P0 = P(1);
    end
    function kppedit(varargin)
        kpp = str2num(varargin{1}.String)*1e-30 / 1e-3;
    end
    function betaedit(varargin)
        gammaK = str2num(varargin{1}.String);
    end
    function zwedit(varargin)
        zw = str2num(varargin{1}.String);
        gfcn = ifftshift(1/zw/sqrt(2*pi)*exp(-1/2*(z2s/zw).^2));
        Fgfcn = fft(gfcn);
        nrm = real(ifft(Fgfcn.*fft(ones(size(z2))))); Fgfcn = Fgfcn./nrm;
        cfcno = cfcn;
        cfcn = @(xi) (ifft(Fgfcn.*fft(xi)));
%         phi = cfcn(phi);
    end
    function hedit(varargin)
        h=str2num(varargin{1}.String);
    end

    function reset(varargin)
        [phia,fa,Pgc,Dv] = Analytic_Calculations();
        F = sqrt(Pgc).*exp(i*phiI); 
    end
    function toggc(varargin)
        gc = 1-gc;
    end
    function phasepert(varargin)
        F = F.*exp(i*0.25*z2/2/Lc);
    end
    function toggletheory(varargin)
        plottheory = 1-plottheory;
    end
    function cplottime(varargin)
        plotinterval=str2num(varargin{1}.String);
    end


    %% Old power function
    function P0 = Steady_State_Power()
        Av = sqrt(Psat)*ones(2*Ndz,1);
        dfigure;
        disptimer=tic;
        maxerr = Inf; kk=1;
        while maxerr>1e-8 || kk<10
            dAvdz     = uddz  (Av,R1,R2);
            dAvdt = c/n*(-dAvdz - aw/2*Av ...
                     + g0/2.*maxLs*( Av -1/Psat*(Av.^2+2*flipud(Av).^2).*Av ));
            maxerr = max(abs(dAvdt*dt./Av)); kk=kk+1;
            Av = Av + dAvdt*dt*0.8;
            if toc(disptimer)>0.5
                plot(abs(Av).^2); drawnow;
                maxerr
                disptimer=tic;
            end
        end        
        P0 = abs(Av).^2;
        close(gcf);
        function dEo = uddz(Ei,R1v,R2v)
            Epi = [Ei(1:Ndz);Ei(Ndz+1)];
            Emi = [Ei(1);flipud(Ei(Ndz+1:end))];
            Epi = [Emi(2)*sqrt(R1v); Epi];
            Emi = [Emi; Epi(end-1)*sqrt(R2v)];
            dEpo = (Epi(2:end)-Epi(1:end-1))/dz;
            dEmo = (Emi(2:end)-Emi(1:end-1))/dz;
            dEo = [dEpo(1:end-1);-flipud(dEmo(2:end))];
        end
    end

    %% New power function (Created by Levi)
    function [Nodes, Power_fct] = SSP(count)
    % Guess Parameters
    xhigh = Lc; xlow = 0;
    fR_const = 1 * Psat;
    
    % Solver
    solinit = bvpinit(linspace(xlow,xhigh,count),@(x) guess(x));
    sol = bvp4c(@(x,y) bvp4ode(x,y), @(ya,yb) bvp4bc(ya,yb), solinit);

    % Point Evaluation
    xint = linspace(xlow, xhigh, Ndz);
    Sxint = deval(sol,xint);

    % Row to Column
    xint = transpose(xint);
    Sxint = transpose(Sxint);
    
    % Reshaping and changing x-axis from distance to number of nodes
    Node_xint = (Ndz/xhigh)*xint;
    Node_Extend_xint = zeros(2*Ndz,1); Node_Extend_xint(1:Ndz,1) = Node_xint; Node_Extend_xint(Ndz+1:2*Ndz,1) = 2000 + Node_xint;
    plt_Array = zeros(2*Ndz,1); plt_Array(1:Ndz,1) = Sxint(:,1); plt_Array(Ndz+1:2*Ndz,1) = flipud(Sxint(:,2));

    Nodes = Node_Extend_xint;
    Power_fct = plt_Array;
        function g = guess(x)
            omega = (1/xhigh)*log(10);
            g = [fR_const*exp(omega*x)
                    0.2*fR_const*exp(omega*x)];
        end
        function dydx = bvp4ode(x,y)
            dydx = [g0*(1-(1/Psat)*(y(1)+2*y(2)))*y(1)-aw*y(1) 
                -(g0*(1-(1/Psat)*(y(2)+2*y(1)))*y(2)-aw*y(2))];
        end
        function res = bvp4bc(ya,yb)
            res = [ya(1)-R1*ya(2) yb(2)-R2*yb(1)];
        end
    end

    %% Analytical Calculations
    function [phia,fa,Pgc,betaeff] = Analytic_Calculations()
        kppeff = kpp - 2*Psat*(T2^2).^gc*gammaK;
        betaeff = kppeff*(c/n)^3;
        gm = g0/(2*Psat)*(c/n)^2*(T1+1/2*T2)*(P(end)-P(1))/(4*Lc*P(1));
%         gm= g0/(2*Psat)*(c/n)^2*(T1+1/2*T2)*1/(4*Lc)*(exp(2*Lc*am)-1);
        r = g0/(2*Psat)*c/n*3*mean(P)/P(1);
        Aa = sqrt(P0)/sqrt(1+gm/2/r);
        
        gammaKp = gammaK*c/n*3*Km;
        Dg=2*g0*T2^2*(c/n)^3*gc;
%         beta = beta - 1*nK*Dg/r/2;
        
        phia = 1/2*gm/betaeff*Aa^2*(z2-Lc).^2;
        fa = -c/(4*pi*n)*gm/(betaeff/2)*abs(Aa).^2.*(z2-Lc);
        fa = fa -1/(24*pi)*gm/(betaeff/2)*Lc^2*gm*abs(Aa).^4;
        fa = fa -1/(2*pi)*gammaK*c/n*3*Km*abs(Aa).^2;
        
        Pgc = P0 + 1/r*(betaeff/2*1/2*gm/betaeff*Aa^2*2 + ...
                       -Dg/4*(1/2*gm/betaeff*Aa^2*2*(z2-Lc)).^2);
                   
        stabnum = P0*Dg/r*(gm*Lc/2/(betaeff/2-gammaKp*Dg/4/r)).^2
        BW = 1/(24*pi)*2*Lc*am*(g0-aw-am)/kppeff*(T1+1/2*T2)
        
%         phia = pi*n/(2*c*Lc)*BW*(z2-Lc).^2;
    end

end

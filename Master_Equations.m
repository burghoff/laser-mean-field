function soln=Master_Equations(p)
% Simulate a QCL (or any other laser) using the master equations used in 
% in Burghoff's "Unraveling the origin of frequency modulated combs using active cavity mean-field theory"
% (Optica 2020, preprint at ArXiv:2006.12397).
% This reduces the Maxwell-Bloch equations of a laser down to two equations on a fine time grid.
% A GUI allows certain parameters to be adjusted on the fly if desired.

% It also produces a theoretical plot for the theoretical form of the extendon
% (equation (7), currently only valid for a Fabry-Perot cavity with either
% R1=1 or R2=1). To save simulation time it initializes the field with this
% as well.
%
% Input: a parameter structure generated using NLFM_Params
% Output: a solution structure that stores the field evolution
% Example 1: s=Master_Equations(NLFM_Params('numTr',5000,'ampphase',1));
%            Runs for 5000 round trips in amplitude-phase mode, which has
%            less numerical diffusion (false gain curvature) but is more
%            prone to divergence.
% Example 2: s=Master_Equations(NLFM_Params('numTr',5000,'ampphase',0,'dtperTr',4000));
%            Runs for 5000 round trips in linear mode, which has a lot of
%            numerical diffusion (false gain curvature). A finer grid of 4000 nodes
%            needs to be used for it to be stable. Note that the mean-field
%            version does not have these issues, as it is split-step based.
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
% Nt = round(numTr/h);                         % Total number of time steps
Nt = round(numTr* 2*Lc/(c/n) /dt);
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
    fh=dfigure('Position',[399 285 1009 437]); ml=-Inf;; ml=-Inf;
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

%% Primary simulation code
[phia,fa,Pa] = Analytic_Calculations();
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
F = sqrt(Pa).*exp(i*phiI);                                % slow-intensity envelope
[Ep0,Em0]=toNormal(F.*sqrt(K));
% Ep0 = Ep0(1:end-1); Em0 = Em0(2:end);
Ep0(end)=Em0(end)/sqrt(R2); Em0(1)=Ep0(1)/sqrt(R1);
Ap0 = abs(Ep0);   Am0= abs(Em0);
pp0 = unwrap(angle(Ep0))/1.0; pm0 = unwrap(angle(Em0))/1.0;

Ap=Ap0; Am=Am0; pp=pp0; pm=pm0;
Ep=Ep0; Em=Em0;

% Things to save
Nout = min(maxsave,Nt); iio=round(linspace(1,Nt,Nout));
[aF,aA,aphi,af,aS] = deal(zeros(length(F),Nout));
[ts,ashft]=deal(zeros(1,Nout));

[Po,po,ts,errs] = deal(zeros(Nt,1));
Po(1)=abs(Ap0(end)).^2; po(1)=angle(pp0(end));
for ii=2:Nt
    Pp = abs(Ep).^2; 
    Pm = abs(Em).^2; %pm = unwrap(angle(Em));
    
    if p.ampphase
        [dApdz,dAmdz]     = uddz  (Ap,Am,R1,R2);
        [dAp2dz2,dAm2dz2] = ud2dz2(Ap,Am,R1,R2);
        [dppdz,dpmdz]     = uddz  (pp,pm,1,1);
        [dpp2dz2,dpm2dz2] = ud2dz2(pp,pm,1,1);

        dApdt0 = -c/n*dApdz;
        dAmdt0 =  c/n*dAmdz;
        dppdt0 = -c/n*dppdz;
        dpmdt0 =  c/n*dpmdz;
        dAp2dt20 = (c/n)^2*dAp2dz2;
        dAm2dt20 = (c/n)^2*dAm2dz2;
        dpp2dt20 = (c/n)^2*dpp2dz2;
        dpm2dt20 = (c/n)^2*dpm2dz2;

        dApdt = c/n*(-dApdz - 1*kpp/2*(2*dApdt0.*dppdt0 + Ap.*dpp2dt20) ...
                     - aw/2*Ap ...
                     + g0/2.*( Ap -1/Psat*(Pp+2*Pm).*Ap - T2*dApdt0 + gc*T2^2*(dAp2dt20-dppdt0.^2.*Ap)+...
                     1*1/Psat.*((3*T1+11/2*T2)*dAmdt0.*Am.*Ap + ...
                                useN2N3*((T1+5/2*T2)*Am.*Am.*dApdt0 + (2*T1+4*T2)*Ap.*Ap.*dApdt0))));
        dAmdt = c/n*(+dAmdz - 1*kpp/2*(2*dAmdt0.*dpmdt0 + Am.*dpm2dt20) ...
                     - aw/2*Am ...
                     + g0/2.*( Am -1/Psat*(Pm+2*Pp).*Am - T2*dAmdt0 + gc*T2^2*(dAm2dt20-dpmdt0.^2.*Am)+...
                     1*1/Psat.*((3*T1+11/2*T2)*dApdt0.*Ap.*Am + ...
                                useN2N3*((T1+5/2*T2)*Ap.*Ap.*dAmdt0 + (2*T1+4*T2)*Am.*Am.*dAmdt0))));
        dppdt = c/n*(-dppdz + kpp/2*(dAp2dt20./Ap - dppdt0.^2) ...
                     -i*gammaK*(Pp+2*Pm) ...
                     +g0/2.*(1*-T2*dppdt0 + gc*T2^2*(2*dApdt0.*dppdt0./Ap + dpp2dt20) + ...
                     1/Psat.*(-(T1+1/2*T2)*Am.*Am.*dpmdt0 + ...
                                useN2N3*((T1+5/2*T2).*Am.*Am.*dppdt0 + T2*Ap.*Ap.*dppdt0))));
        dpmdt = c/n*(+dpmdz + kpp/2*(dAm2dt20./Am - dpmdt0.^2) ...
                     -i*gammaK*(Pm+2*Pp) ...
                     +g0/2.*(1*-T2*dpmdt0 + gc*T2^2*(2*dAmdt0.*dpmdt0./Am + dpm2dt20) + ...
                     1/Psat.*(-(T1+1/2*T2)*Ap.*Ap.*dppdt0 + ...
                                useN2N3*((T1+5/2*T2).*Ap.*Ap.*dpmdt0 + T2*Am.*Am.*dpmdt0))));         

        Ap = Ap + dApdt*dt;
        Am = Am + dAmdt*dt;
        pp = pp + dppdt*dt;
        pm = pm + dpmdt*dt;

        Ep = Ap.*exp(i*pp); Em = Am.*exp(i*pm);
    else
        [dEpdz,dEmdz]     = uddz(Ep,Em,R1,R2);
        [dEp2dz2,dEm2dz2] = ud2dz2(Ep,Em,R1,R2);

        dEpdt0 = -c/n*dEpdz;
        dEmdt0 =  c/n*dEmdz;
        dEpdt = c/n*(-dEpdz + i*kpp/2*(c/n)^2*dEp2dz2 - aw/2*Ep + ...
                         -i*gammaK*(Pp+2*Pm).*Ep + ...
                         g0/2.*( Ep -1/Psat*(Pp+2*Pm).*Ep - T2*dEpdt0 + gc*T2^2*(c/n)^2*dEp2dz2 + ...
                     1/Psat.*((2*T1+3*T2)*conj(dEmdt0).*Em.*Ep + (T1+5/2*T2)*conj(Em).*dEmdt0.*Ep + ...
                              useN2N3*(T1+5/2*T2)*conj(Em).*Em.*dEpdt0 + useN2N3*(T1+5/2*T2)*conj(Ep).*dEpdt0.*Ep + ...
                              useN2N3*(T1+3/2*T2)*Ep.*conj(dEpdt0).*Ep)));
        dEmdt = c/n*(+dEmdz + i*kpp/2*(c/n)^2*dEm2dz2 - aw/2*Em + ...
                     -i*gammaK*(Pm+2*Pp).*Em + ...
                         g0/2.*( Em -1/Psat*(Pm+2*Pp).*Em - T2*dEmdt0 + gc*T2^2*(c/n)^2*dEm2dz2 + ...
                     1/Psat.*((2*T1+3*T2)*conj(dEpdt0).*Ep.*Em + (T1+5/2*T2)*conj(Ep).*dEpdt0.*Em + ...
                              useN2N3*(T1+5/2*T2)*conj(Ep).*Ep.*dEmdt0 + useN2N3*(T1+5/2*T2)*conj(Em).*dEmdt0.*Em + ...
                              useN2N3*(T1+3/2*T2)*Em.*conj(dEmdt0).*Em)));

        Ep = Ep + dEpdt*dt;
        Em = Em + dEmdt*dt;
    end
    
%     Eps = [Eps(:,end),Ep];
%     Ems = [Ems(:,end),Em];

    Po(ii)=abs(Ep(end)).^2;
    tmp=unwrap([po(ii-1),angle(Ep(end))]);
    po(ii)=tmp(2);
    ts(ii)=ts(ii-1)+dt;
    
%     Pp = 
  
    if plotprogress & toc(disptimer)>2
        pr=round(linspace(1,ii,min(3000,ii)))'; pr=pr(1:end-1);
        makeplot(pr)
        if pushed==1, break; end
    end
    if any(iio==ii)
        myi = find(iio==ii);
        phi = unwrap(angle(toFlipped(Ep,Em)));
%         dphiend = round((phi(end)-phi(1))/(2*pi))*2*pi;
%         phi = phi - dphiend/(2*Lc)*z2;
%         [~,l6] = min(abs(z2-0.75*2*Lc));
%         [~,ml]=max(phi/sign(kpp));
        aE   (:,myi)=toFlipped(Ep,Em);
        aphi (:,myi)=phi;
        aA   (:,myi)=abs(toFlipped(Ep,Em));
        af   (:,myi)=diff([phi(end);phi])/dz*(-c/n)/(2*pi)/1e12;
        aS   (:,myi)=fftshift(fft(toFlipped(Ep,Em)));
%         ashft(myi)=l6-ml;
    end
end

% Construct solution
soln=struct();
sr = [1:myi];  % in case we ended early
soln.z=z2;                                                          % position
soln.P=P;                                                           % steady state power function
soln.t=ts(iio(sr));                                                 % times (s)
soln.E=aE(:,sr);                                                    % envelopes
soln.phi=aphi(:,sr);                                                % inst. phases
soln.A=aA(:,sr);                                                    % amplitudes
soln.f=af(:,sr);                                                    % inst. frequency
% soln.shft=ashft(sr);                                                % shifts to move the jump to 75%
soln.S=aS(:,sr);                                                    % spectrum
soln.Sf=-c/n*fss;                                                   % spectrum frequencies
% soln.Ls=real(Ls(-c/n*fss*2*pi));                                    % gain lineshape

% Last 10% metrics 
% [soln.pn10] = deal(zeros(size(fss)));                             
% for jj=1:length(soln.Sf)
%     phiS = unwrap(angle(soln.S(jj,:)));
%     gr = [round(myi*9/10):myi];
%     dphi = phiS-polyval(polyfit(soln.t(gr),phiS(gr),1),soln.t);
%     soln.pn10(jj)=mean(abs(dphi(gr)).^2);                           % phase noise
% end
% soln.g10  = real(sqrt(1-soln.pn10));                                % g^(1) coherence (spectrally-resolved)
% soln.aPS10 = mean(abs(soln.S(:,gr)).^2,2);                          % average power spectrum
% soln.ag10 = sqrt(sum(soln.g10.^2.*soln.aPS10)./sum(soln.aPS10));    % average g^(1), weighted by power
                 

soln.zn = z;                                                        % normal cavity position (m)
soln.toNormal = @toNormal;                                          % function that returns forward and backward waves          
soln.params=p;                                                      % input parameters
function [Epi,Emi] = toNormal(Ei)
    Epi = [Ei(1:Ndz);NaN*Ei(Ndz+1)];
    Emi = [Ei(1)*NaN;flipud(Ei(Ndz+1:end))];
end
function Eo = toFlipped(Epi,Emi)
    Eo = [Epi(1:end-1);flipud(Emi(2:end))];
end

[phia,fa,Pa,Deff,fceo_XS,fceo_Kerr] = Analytic_Calculations();
soln.analytic = struct('phi',phia,'f',fa,'P',Pa,'betaeff',Deff,'fceo_XS',fceo_XS,'fceo_Kerr',fceo_Kerr);

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
        [phia,fa,Pa,Dv] = Analytic_Calculations();
        F = sqrt(Pa).*exp(i*phiI); 
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
    function [phia,fa,Pa,betaeff,fceo_XS,fceo_Kerr] = Analytic_Calculations()
        kppeff = kpp - 2*Psat*(T2^2).^gc*gammaK;
        betaeff = kppeff*(c/n)^3;
        gm = g0/(2*Psat)*(c/n)^2*(T1+1/2*T2)*(P(end)-P(1))/(4*Lc*P(1));
%         gm= g0/(2*Psat)*(c/n)^2*(T1+1/2*T2)*1/(4*Lc)*(exp(2*Lc*am)-1);
        r = g0/(2*Psat)*c/n*3*mean(P)/P(1);
        Aa = sqrt(P0)/sqrt(1+gm/2/r);
        
        gammaKp = gammaK*c/n*3*Km;
        Dg=2*g0*T2^2*(c/n)^3*gc;
        
        phia = 1/2*gm/betaeff*Aa^2*(z2-Lc).^2;
        fa = -c/(4*pi*n)*gm/(betaeff/2)*abs(Aa).^2.*(z2-Lc);
        fceo_XS = -1/(24*pi)*gm/(betaeff/2)*Lc^2*gm*abs(Aa).^4;
        fceo_Kerr = -1/(2*pi)*gammaK*c/n*3*Km*abs(Aa).^2;
        fa = fa + fceo_XS + fceo_Kerr;
        
%         Pgc = P0 + 1/r*(betaeff/2*1/2*gm/betaeff*Aa^2*2 + ...
%                        -Dg/4*(1/2*gm/betaeff*Aa^2*2*(z2-Lc)).^2);
        Pa = Aa.^2 - Dg/(4*r)*(Aa^2*gm*(z2-Lc)/betaeff).^2;       % bug fixed 1.20.2021
                   
        stabnum = P0*Dg/r*(gm*Lc/2/(betaeff/2-gammaKp*Dg/4/r)).^2
        BW = 1/(24*pi)*2*Lc*am*(g0-aw-am)/kppeff*(T1+1/2*T2)
    end

    function [dEpo,dEmo] = uddz(Epi,Emi,R1v,R2v)
        Epi = [Emi(2)*sqrt(R1v); Epi];
        Emi = [Emi; Epi(end-1)*sqrt(R2v)];
        dEpo = (Epi(2:end)-Epi(1:end-1))/dz;
        dEmo = (Emi(2:end)-Emi(1:end-1))/dz;
    end
    function [dEpo,dEmo] = ud2dz2(Epi,Emi,R1v,R2v)
        Epi = [Emi(3)*sqrt(R1v);Emi(2)*sqrt(R1v); Epi];
        Emi = [Emi; Epi(end-1)*sqrt(R2v); Epi(end-2)*sqrt(R2v)];
        dEpo = (1*Epi(3:end)-2*Epi(2:end-1)+1*Epi(1:end-2))/dz^2;
        dEmo = (1*Emi(3:end)-2*Emi(2:end-1)+1*Emi(1:end-2))/dz^2;
    end

    function makeplot(pr)
        disptimer = tic;
        ctime = round(ii*dt/Tr)
%         dta/dt
        subplot(2,3,1);
        plot(z/1e-3,abs(Ep).^2,z/1e-3,abs(Em).^2);
%         plot(z/1e-3,Pp,z/1e-3,Pm,z/1e-3,Pp+Pm);
%         plot(z/1e-3,-g0/2/Psat*(T1+1/2*T2)*Pm.*c./n.*dpmdt0,z/1e-3,-g0/2/Psat*(T1+1/2*T2)*Pp.*-c./n.*dppdt0);
%         ml=get(gca,'YLim'); ylim([0,max(ml)]);
        xyt('Position (mm)','Power','');
        subplot(2,3,4);
        pplot = unwrap(angle(toFlipped(Ep,Em)));
        [ppplot,pmplot] = toNormal(pplot-mean(pplot));
        pplot0 = unwrap(toFlipped(pp0,pm0));
        [ppplot0,pmplot0] = toNormal(pplot0-mean(pplot0));
%         [dppdz,dpmdz]     = uddz(ppplot,pm2,R1,R2); dppdz(1)=NaN; dpmdz(end)=NaN;
        plot(z/1e-3,ppplot,z/1e-3,pmplot); 
        hold all; plot(z/1e-3,ppplot0,z/1e-3,pmplot0); hold off;
        xyt('Position (mm)','Phase','');
        
        a2=subplot(2,3,5);
        ti = dt*(pr-1);
        fi = (po(pr+1)-po(pr))/dt/2/pi;
        plot(ti/Tr,fi);  xyt('Time/T_r','\Delta f (Hz)',''); hold all;
        fo2 = reshape([(po(2:ii)-po(1:ii-1))/2/pi/dt;NaN*ones(ceil(ii/round(Tr/dt))*round(Tr/dt)-ii+1,1)],...
                      round(Tr/dt),ceil(ii/round(Tr/dt)));
        plot(dt*size(fo2,1)*[0:size(fo2,2)-1]/Tr,max(fo2,[],1,'includenan'),...
             dt*size(fo2,1)*[0:size(fo2,2)-1]/Tr,min(fo2,[],1,'includenan')); hold off;          
                  
        a5=subplot(2,3,2);
        plot(ti/Tr,Po(pr)); xyt('Time/T_r','Power','');
        hold all;
        Po2 = reshape([Po(1:ii);NaN*ones(ceil(ii/round(Tr/dt))*round(Tr/dt)-ii,1)],...
                      round(Tr/dt),ceil(ii/round(Tr/dt)));
        plot(dt*size(Po2,1)*[0:size(Po2,2)-1]/Tr,max(Po2,[],1,'includenan'),...
             dt*size(Po2,1)*[0:size(Po2,2)-1]/Tr,min(Po2,[],1,'includenan')); hold off;
%         linkaxes([a2,a5],'x');
        
        a3=subplot(2,3,6);
        pr=max(1,[ii-round(5*Tr/dt):ii])'; pr=pr(1:end-1);
        ti = dt*(pr-1);
        fi = (po(pr+1)-po(pr))/dt/2/pi;
        plot(ti/Tr,fi); xyt('Time/T_r','\Delta f (Hz)','');
        a6=subplot(2,3,3);
        plot(ti/Tr,Po(pr)); xyt('Time/T_r','Power','');
        drawnow;
%         linkaxes([a3,a6],'x');
    end

end

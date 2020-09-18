function os = NLFM_Params(varargin)
% Generates an input structure for a nonlinear FM simulation

% Device params
os.J = 1100;                        % current density, A/cm^2
os.kpp = -2000;                      % dispersion, fs^2/mm
os.gammaK = 0;                      % Kerr nonlinearity, m/V^2
os.R1=0.09;                         % mirror 1 reflectivity
os.R2=1;                            % mirror 2 reflectivity
os.lm0=8e-6;                        % Wavelength, m
os.Lc=4e-3;                         % Laser cavity length, m
os.df=0e12;                         % Frequency detuning, Hz

% QCL parameters
os.mu=2.3e-9;                       % Dipole matrix element, nm
os.T1=0.4e-12;                      % Population lifetime, s
os.T2=50e-15;                       % Coherence lifetime, s
os.n=3.3;                           % Group refractive index
os.Lmod=580e-10;                    % Module length, m

% Waveguide parameters
os.aw=4;                            % Waveguide loss (cm^-1)
os.Gamma=1;                         % Gamma (waveguide overlap)

% Simulation params
os.h = 1;                           % split-step timestep (fraction of round trip time)
os.dtperTr = 1000;                  % nodes per round trip
os.numTr = 10000;                   % number of round trip times to simulate
os.Crnt=1;                          % Courant number
os.gc=1;                            % gain curvature enabled?
os.Ls='lorentzian';                 % lineshape function ('lorentzian' or 'parabolic')
os.useN2N3=0;                       % use all the nonlinear terms? (requires smaller h, 0.125 sufficient for default settings)
os.initphi='fundamental';           % phase initialization ('fundamental','harmonicN','cosineN','rand')
os.maxsave=10000;                   % number of states to save (1 gives final result only)
os.useLinPwrApprox = 0;             % approximate cross-steepening terms using the linear power approximation? (leads to the NLSE)

% Plot parameters
os.plotprogress=1;                  % plot progress?
os.plottheory = 1;                  % plot theoretical result?
os.plotinterval = 0.5;                  % how frequently to plot

fields = fieldnames(os);
for ii=1:length(varargin)
    for jj=1:length(fields)
        if ischar(varargin{ii}) & isequal(varargin{ii},fields{jj})
            os.(fields{jj}) = varargin{ii+1};
        end
    end     
end

end
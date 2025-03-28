clear all
close all
clc

addpath('Functions');


%% PHYSICAL CONSTANTS
global c; c = 2.99792458e8; % Speed of radiowaves in vacuum expressed in [m.s^-1]


%% SPECIFY THE RADAR PARAMETERS
%  Waveform parameters
f0          =   10.00e9                            ; % Operating frequency expressed in [Hz]
B           =  300.00e6                            ; % Transmitted bandwidth expressed in [Hz]
k0          = 2*pi*f0/c                            ; % Wavenumber in the vacuum at f0 expressed in [m^-1]

%  LWA parameters
epsilon_r   =   10.00                              ; % Permittivity of the dielectric material filling up the waveguide expressed in []
h           =    c/(2*(f0-3*B/4) * sqrt(epsilon_r)); % Height of the waveguide expressed in [m]
P           =    0.40  * c/f0                      ; % Slots spacing expressed in [m]
N_slots     =   32                                 ; % Number of slots
alpha       =    0.015 *   k0                      ; % Leakage rate expressed in [rad.m^-1]

%  Notes:
%  - The cut-off frequency (fc) is obtained when beta = sqrt((2*pi*fc/c)^2 * epsilon_r - (pi/h)^2) = 0, i.e. beta becomes imaginary below fc.
%    Therefore, h = c/(2*fc*sqrt(epsilon_r)).
%
%  - The efficiency is defined as 1 - exp(-2*alpha*N_slots*P) and is generally chosen between 0.80 and 0.95.
%
%  - The parameter alpha can be adapted by changing the dimensions of the slots.


%% ANTENNA PATTERN vs. ANGLE
%  Initializations
Legend  =                             [];
f       = linspace(f0-B/2, f0+B/2,    9); % Frequency expressed in [Hz]
u       = linspace( -1.00,  +1.00, 1001); % Direction cosine expressed in []
WA      =  complex(zeros(1,  length(u))); % LWA transfer function

%  Create the figure
Figure  = figure(1); set(Figure, 'Position', [100 100 400 350]); hold on; grid on;
xlabel('asin$(u)$ [deg]'                   , 'Interpreter', 'Latex');
ylabel('$|\bf{W}^T(f)\bf{A}(f, u)|$ [dB]'  , 'Interpreter', 'Latex');
 title('LWA diagram'                       , 'Interpreter', 'Latex');

% Compute and display the antenna pattern as function of the angle and parametrized by the frequency
for kf = 1:length(f)
    
    %  LWA model
    for ku = 1:length(u), [WA(ku)] = LWA(u(ku), f(kf), h, P, N_slots, alpha, epsilon_r);
    end
        
    %  Display the antenna diagram
    plot(asind(u), 20*log10(abs(WA)));
    
    %  Update the legend
    Legend = [Legend ; sprintf('$f = %5.2f$~GHz', f(kf)/1e9)]; %#ok<AGROW>
        
end

%  Create the figure
legend(Legend, 'Interpreter', 'Latex', 'Location', 'East'); axis([0 90 -15.00 +15.00]);


%% OBSERVATION OF THE RECEIVED SIGNAL
%  Specify the target
u	= 0.70; % Direction cosine of the target expressed in []

%  Specify the simulation parameters
Fs	= f0+B; % Sampling frequency expressed in [Hz]
N 	= 1024; % Number of samples

%  Define the time vector expressed in [s]
t   = (0:N-1)/Fs;

%  Specify the transmitted signal
s   = exp(1j * 2*pi * ((f0 - B/2) * t + B * t.^2/(2* t(end)))); % Chirp

%  Compute the LWA response to the target backscattered signal
f   =                    (0:N-1)/N * Fs  ;
WA  =  complex(zeros(1      , length(f))); % LWA transfer function
W   =  complex(zeros(N_slots, length(f))); % Waveguide transfer function
 A  =  complex(zeros(N_slots, length(f))); % Slotted array response in the direction cosine u

for kf = 1:length(f), [WA(kf), W(:, kf), A(:, kf), p] = LWA(u, f(kf), h, P, N_slots, alpha, epsilon_r);
end

%  Received signal in a noiseless case
r = ifft(fft(s).* WA.^2);

%  Display the power spectrums
Figure  = figure(2); set(Figure, 'Position', [100 100 400 350]); set(gca, 'TickLabelInterpreter', 'Latex'); hold on; grid on;
xlabel('$f$ [GHz]',                       'Interpreter', 'Latex');
ylabel('psd [dB ]',                       'Interpreter', 'Latex');
 title('Transmitted and received signal', 'Interpreter', 'Latex');
 
plot((0:N-1)/N * Fs /1e9, 20*log10(abs(fft(s))), 'b');
plot((0:N-1)/N * Fs /1e9, 20*log10(abs(fft(r))), 'r');
plot((0:N-1)/N * Fs /1e9, 20*log10(abs(   WA )), 'k');

legend(        'Transmitted signal                             ',               ...
       sprintf('Received signal ($u = %.2f$, namely $%.1f$ deg)', u, asind(u)), ...
       sprintf('LWA response   ~($u = %.2f$, namely $%.1f$ deg)', u, asind(u)), 'Interpreter', 'Latex');

axis([(f0-B)/1e9 (f0+B)/1e9 -20 inf]);

%  Display the time sequences
Figure  = figure(3); set(Figure, 'Position', [100 100 400 350]); set(gca, 'TickLabelInterpreter', 'Latex'); hold on; grid on;
xlabel('$t$ [s]',                         'Interpreter', 'Latex');
ylabel('Amplitude',                       'Interpreter', 'Latex');
 title('Transmitted and received signal', 'Interpreter', 'Latex');
 
plot(t, real(s        ), 'b');
plot(t, real(r/N_slots), 'r');

legend('$\\Re(s)$', '$\\Re(r)$', 'Interpreter', 'Latex');


%% COMPUTE THE ULA PARAMETERS
%  Compute the radiated power by each channel of the ULA is therefore Pt
Pt = 1 / N_slots;

%  Compute the average transmit gain of the ULA in the direction of the target
sqrt_g = 0;

for v = linspace(sind(30.00), sind(60), length(t))
    
    [~, ~, Av] = LWA(v, f0, h, P, N_slots, alpha, epsilon_r); % Beamforming coefficients
    [~, ~, Au] = LWA(v, f0, h, P, N_slots, alpha, epsilon_r); % ULA response to the target
    
    sqrt_g     = sqrt_g  + abs(Av' * Au);
    
end

Gt = sqrt_g^2 / length(v);


%% COMPARE THE CRAMER-RAO BOUND FOR ANGLE ESTIMATION WITH A LWA AND WITH AN ULA
%  Specify the simulation parameters
M               = 10000                                        ; % Number of Monte-Carlo runs
sigma_squared   = 10.^(linspace(70.00, 20.00,  51)/10)         ; % Noise power
SNR             = 10*log10(Pt * Gt * length(v)./ sigma_squared); % Equivalent SNR for a ULA-based radar expressed in [dB]
Grid            =      linspace(sind(20.00), sind(60.00), 1001); % Grid containing the direction cosine hypothesises expressed in []

%  Declaration of variables
% MSE_LWA = zeros(size(SNR));
CRB_LWA = zeros(size(SNR));
CRB_ULA = zeros(size(SNR));

%  Loop over SNR values
for k = 1:length(SNR)
	
    k

    CRB_LWA(k) = Compute_CRB_LWA(s, Fs, sigma_squared(k), p, WA, A, W);
	CRB_ULA(k) = Compute_CRB_ULA(u, P, N_slots, SNR(k)               );
    
    %  Monte-Carlo simulations
    u_hat = zeros(1, M); 
    for m =       1: M
        
        %  Generate the noise vector
        epsilon = sqrt(sigma_squared(k)/2) * (randn(size(s)) + 1j * randn(size(s)));
        
        % MLE
        [u_hat(m)] = MLE(r + epsilon, s, f, h, P, N_slots, alpha, epsilon_r, Grid);

    end

    %  Compute the MSE
    MSE_LWA(k) = mean((u - u_hat).^2);
    
end


%  Display the performance
Figure  = figure(4); set(Figure, 'Position', [100 100 400 350]); set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'Latex'); hold on; grid on;

plot(SNR           - 20*log10(N_slots), MSE_LWA        , '-b ', 'LineWidth', 1.00, 'MarkerSize', 6.00);
plot(SNR(1:6:end) - 20*log10(N_slots), CRB_LWA(1:6:end), '-m^', 'LineWidth', 1.00, 'MarkerSize', 6.00);
plot(SNR(4:6:end) - 20*log10(N_slots), CRB_ULA(4:6:end), '-ko', 'LineWidth', 1.00, 'MarkerSize', 6.00);

xlabel('$N\,\sigma^2/2\,G_t\,T\,|b|^2$ [dB]', 'Interpreter', 'Latex');
ylabel('E$[(u-\hat{u})^2]$',                  'Interpreter', 'Latex');
legend('MSE LWA', 'CRB LWA', 'CRB ULA',       'Interpreter', 'Latex');
axis([-20 +20 1e-7 1e-1]);

save('Workspace');
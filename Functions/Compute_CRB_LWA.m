function [CRB] = Compute_CRB_LWA(s, Fs, sigma_squared, p, GA, A, G)

% Inputs:
% - s            : Transmitted signal.
% - Fs           : Sampling frequency expressed in [Hz].
% - sigma_squared: Noise power (linear scale).
% - p            : Slots spacing expressed in [m].
%  - GA          : LWA reponse.
%  - G           : Waveguide response.
%  - A           : Slotted array response.
%
% Output:
% - CRB          : CRB on the estimation of u expressed in [].

global c;

%% COMPUTE THE PARTIAL DERIVATIVES
%  Specify the frequency vector
N =      length(s); % Number of samples
f = (0:N-1) * Fs/N; % Frequency vector expressed in [Hz]

%  Specify the wavenumber
k = 2*pi * f / c;

%  Derivative with respect to the direction cosine (u)
dA_du       = 1j * k.* p.* A;  
fft_ds_du   = 2*fft(s).* GA.* sum(G.* dA_du, 1); % Chain rule
ds_du       = ifft(fft_ds_du); % Convert back to time domain


%% COMPUTE THE FISHER INFORMATION AND THEN THE CRB
FIM = 2 * sum(abs(ds_du).^2) / sigma_squared;
CRB = 1 / FIM                               ;
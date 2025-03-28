function [WA, W, A, p] = LWA(u, f, h, P, N, alpha, epsilon_r)

%  Inputs:
%  - u       : Direction cosine expressed in [].
%  - f       : Frequency expressed in [Hz].
%  - h       : Height of the waveguide expressed in [m].
%  - P       : Slots spacing expressed in [m].
%  - N       : Number of slots.
% - alpha    : Leakage rate expressed in [rad.m^-1].
% - epsilon_r: Permittivity of the dielectric material filling up the waveguide expressed in [].
%
%  Outputs:
%  - W       : Waveguide response.
%  - A       : Slotted array response in the direction cosine u.
%  - WA      : LWA reponse in the direction cosine u.
%  - p       : Positions of the slots along the waveguide shallowness expressed in [m.

global c;

%% COMPUTE ADDITIONNAL ANTENNA PARAMETERS
k0      = 2*pi*f/c                         ; % Wavenumber in the vacuum expressed in [m^-1]
beta    = sqrt(k0^2 * epsilon_r - (pi/h)^2); % Phase constant expressed in [rad.m^-1]
p       = (0:N-1).' * P                    ; % Slots positions expressed in [m]


%% LWA MODEL
%  Waveguide response
W = exp( -alpha * p - 1j * beta * p) * isreal(beta);

if norm(W) ~= 0, W = W/norm(W);
end

% Slotted array response in the direction cosine u
A = exp(1j * k0*u*p);                                   

% LWA reponse in the direction cosine u
WA = W.' * A;
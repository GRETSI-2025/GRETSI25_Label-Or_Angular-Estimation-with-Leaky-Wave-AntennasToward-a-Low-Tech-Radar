function [u_hat] = MLE(r, s, f, h, P, N, alpha, epsilon_r, Grid)

%  Inputs:
%  - r        : Received signal.
%  - s        : Transmitted signal.
%  - f        : Transmitted frequency expressed in [Hz].
%  - h        : Height of the waveguide expressed in [m].
%  - P        : Slots spacing expressed in [m].
%  - N        : Number of slots.
%  - alpha    : Leakage rate expressed in [rad.m^-1].
%  - epsilon_r: Permittivity of the dielectric material filling up the waveguide expressed in [].
%  - Grid     : Grid containing the direction cosine hypothesises expressed in [].
%
%  Outputs:
%  - u_hat    : Estimate of the direction cosine expressed in [].

persistent r_hat;

%% CONSTRUCT THE REPLICAE OF THE RECEIVED SIGNAL
if isempty(r_hat)
    
    %  Initialization
    r_hat = complex(zeros(length(Grid), length(f)));
    
    %  Loop over the grid
    for g = 1:length(Grid)
    
        %  Compute the LWA response to the target backscattered signal
        GA  =  complex(zeros(1,  length(f))); % LWA transfer function

        for kf = 1:length(f), [GA(kf)] = LWA(Grid(g), f(kf), h, P, N, alpha, epsilon_r);
        end

        %  Construct the replicae of the received signal
        r_hat(g, :) = ifft(fft(s).* GA.^2);
    
    end
end


%% MAXIMUM LIKELYHOOD ESTIMATOR
%  Initialization
log_likelihood = zeros(size(Grid));

%  Loop over the grid
for g = 1:length(Grid)  
    
    %  Compute the log-likelyhood function
    log_likelihood(g) = -norm(r - r_hat(g, :))^2;

end

% Maximize the log-likelyhood
[~, g_hat] =  max(log_likelihood);
    u_hat  = Grid(g_hat         );
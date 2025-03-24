function [CRB] = Compute_CRB_ULA(u, P, N, SNR)

% Inputs:
% - u  : Direction cosine expressed in [].
% - P  : Channels spacing expressed in [m].
% - N  : Number of channels.
% - SNR: Signal-to-Noise Ratio expressed in [dB].
%
% Output:
% - CRB: CRB on the estimation of u expressed in [].
%
%  The equations are detailed in:
%  Document [1]: Overview of Generalized Monopulse Estimation - Ulrich Nickel, IEEE A&E Systems Magazine, Volume 21, (6), June 2006
%  Document [2]: Cramer-Rao Bounds for Estimating Range, Velocity, and Direction with a Sensor Array - Aleksandar Dogandzic, and Arye Nehorai, Proceedings of the 2000 IEEE Sensor Array and Multichannel Signal Processing Workshop (March 2000)
%  Document [3]: Derivation of Monopulse Angle Accuracy for Phased Array RADAR to achieve Cramer-Rao lower Bound - Ryuhei Takahashi, Kazufumi Hirata, Teruyuki Hara, and Atsushi Okamura, ICASSP (25-30 March 2012)


%% SPECIFY THE PARAMETERS
%  Specify the channels positions expressed in [m]
p = (1:N).' * P;
p = p - mean(p);

% Specify the noise and target power
Q         = eye(N     );
b_squared = 10^(SNR/10);


%% COMPUTE THE CRB
%  CRB in the general case (cf. Document [1] - Equation (19) or Document [2] - Equation (3.1a))
a 	= exp(1j * 2*pi * p*u)  ;
a_u	=     1j * 2*pi * p .* a;

CRB = 1/(2 * b_squared) * inv(real(a_u' * inv(Q) * a_u - (a_u' * inv(Q) * a) * (a' * inv(Q) * a_u) / (a' * inv(Q) * a))); %#ok<MINV

%  Note:
%  The CRLB can also be written when Q = Id (cf. Document [3] - Equation (9))
%  CRB = 1/(2*pi)^2 * 1/(sum((p - mean(p)).* (p - mean(p)))) * (1/2) * 1/b_squared
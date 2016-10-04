% Copyright (c) 2011, Ilker R Capoglu 
% Calculates the center frequency (f0), center wavelength, 
% time constant (tau), and -40dB wavelength of the 
% modulated Gaussian waveform that has the specified cutoff wavelengths.

attn_at_cutoff_lambdas = -40;  % attenuation at the cutoff frequencies, in dB
min_lambda_cutoff = 400e-9;  % minimum wavelength in free space (in m)
max_lambda_cutoff = 700e-9;  % maximum wavelength in free space (in m)

c=299792458; %speed of light in free space (in m)
min_k = 2*pi/max_lambda_cutoff;
max_k = 2*pi/min_lambda_cutoff;
min_w = min_k*c;
max_w = max_k*c;
tau_factor = sqrt(-log(10^(attn_at_cutoff_lambdas/20))*2); %number of 1/tau's to reach the cutoff

% Modulated Gaussian pulse: -20dB at w=+-2.15/tau, -40dB at w=+-3.035/tau
w0 = double((min_w+max_w)/2);
half_BW = (max_w-min_w)/2;
tau = tau_factor/half_BW;

dispformat = '%11.5e';
disp(['Center frequency (Hz): ' num2str(w0/2/pi,dispformat)]);
disp(['Center wavelength (m): ' num2str(c/(w0/2/pi),dispformat)]);
disp(['Time constant (sec): ' num2str(tau,dispformat)]);

tau_factor_40dB = sqrt(-log(0.01)*2); %number of 1/tau's to reach -40dB
max_40dB_freq = (w0+tau_factor_40dB*1/tau)/2/pi;
min_40dB_lambda = c/max_40dB_freq;
disp(['[Min. -40dB wavelength (m): ' num2str(min_40dB_lambda) ']']);
disp(['[Recommended spatial step size (m): ' num2str(min_40dB_lambda/15) ']']);

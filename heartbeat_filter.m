function Hd = chebyshev_IIR2
%CHEBYSHEV_IIR
% Chebyshev Type I Bandpass filter designed using FDESIGN.BANDPASS.
% All frequency values are in Hz.
Fs = 20;  % Sampling Frequency

Fstop1 = 0.7;        % First Stopband Frequency
Fpass1 = 0.8;         % First Passband Frequency
Fpass2 = 1.9;         % Second Passband Frequency
Fstop2 = 1.95;        % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 3;           % Passband Ripple (dB)
Astop2 = 60;          % Second Stopband Attenuation (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd = design(h, 'cheby1', 'MatchExactly', match);

% [EOF]

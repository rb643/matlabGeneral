function [p, f] = powerSpectrum(TC, samplingRate)

% subtract mean (DC component)

%TC = TC(1:500);
TC = TC(~isnan(TC));
TC = TC - mean(TC);
fs = samplingRate; %% sampling rate
m = length(TC);          % Window length
n = pow2(nextpow2(m));  % Transform length
y = fft(TC,n);           % DFT
f = (0:n-1)*(fs/n);     % Frequency range

% convert to rms units
p = sqrt(2).*abs(y)./length(y);

end



% range = 2:floor(length(p)/2);
% figure;
% plot(f(range),p(range))
% xlabel('Frequency (Hz)')
% ylabel('Power (rms)')
% title('{\bf powerspectrum}')

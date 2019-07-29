% See: http://biomedicalsignalandimage.blogspot.com/2016/02/matlab-code-to-estimate-power-spectrum.html

% DELTA (0.5 to 4 Hz) THETA(4 to 8 Hz), APLA( 8 to 12 Hz),BETA( 12 to 30 Hz),GAMMA( >30 Hz) 

%close all; clear all; clc;
% fs Sampling frequency, positive scalar. Sampling frequency, specified as
% a positive scalar. The sampling frequency is the number of samples per
% unit time. If the unit of time is seconds, the sampling frequency has units of hertz.
fs = EEG.srate; %  3000; 
% sampling rate or frequency
T = 1/fs;
% contains eeg1 and fs
%load('J:\BIOM_Signal_processing\Hw9\EEGsignal_1') 
eeg1 = EEG.data(1,:);

N =length(eeg1);
% find the length of the data per second
ls = size(eeg1); 
% Make time axis for EEG signal
tx =[0:length(eeg1)-1]/fs;

%EEG waveform
fignum = figure; 
subplot (211), 
plot(tx,eeg1); 
xlabel('Time (s)'), 
ylabel('Amplitude (mV)'), 
title('Original EEG signal'); 

% Used to zoom in on single ECG waveform
subplot(212), plot(tx,eeg1);
%fignum = fignum + 1; 
%figure(fignum);
figure;
xlabel('Time (s)'), 
ylabel('Amplitude (mV)'), 
title('Zoom into original EEG signal at 1 to 2 seconds'), 
xlim([1,2]) 
figure;
%number points. Number of DFT points, specified as a positive integer. For
%a real-valued input signal, x, the PSD estimate, pxx has length (nfft/2 +
%1) if nfft is even, and (nfft + 1)/2 if nfft is odd. For a complex-valued
%input signal,x, the PSD estimate always has length nfft. If nfft is
%specified as empty, the default nfft is used.  
NFFT = 2 ^ nextpow2(N);  
Y = fft(eeg1, NFFT)/N; %% fft of the EEG signals
f = (fs/2 * linspace(0,1,NFFT/2+1))'; % Vector containing frequencies in Hz
amp =( 2 * abs(Y(1: NFFT/2+1))); % Vector containing corresponding amplitudes
subplot(2,1,1), 
plot (f, amp);
title ('plot single-sided amplitude spectrume of the EEG signal'); 
xlabel ('frequency (Hz)'); ylabel ('|y(f)|');grid on;
%Estimate the power spectrum of the 10-s epoch by computing the periodogram
%the plot using periodogram with no outputs.
periodogram_1 = periodogram(eeg1);
subplot(2,1,2), 
plot (f,periodogram_1);
title('Periodogram power spectral Densisty Estimate of the EEG signal'); 
xlabel ('frequency (Hz)');
ylabel ('Power/Frequency(dB)');
grid on;
%%% plot the periodogram with the length of EEG signal
nfft= length(eeg1);
periodogram_2 =periodogram(eeg1,[],nfft);
f1 = (fs/2 * linspace(0,1,nfft/2+1))'; % Vector containing frequencies in Hz
figure;
subplot(2,2,1), 
plot(f1,periodogram_2);
title('Periodogram power spectral Densisty Estimate of the EEG signal and length of EGG signal'); 
xlabel('frequency (Hz)');ylabel('Power/Frequency(dB)');
grid on;%%
%The signal is 30001 samples in length. Obtain the periodogram using the
%default rectangular window and DFT length. The DFT length is the next power of two greater than the signal length, or 32768 points. 
%Because the signal is real-valued and has even length, the periodogram is one-sided and there are 512/2+1 points.
[pxx,w] = periodogram(eeg1);
subplot(2,2,2), 
plot(w,10*log10(pxx)); 
title ('Periodogram using the default rectangular window and DFT length'); 
xlabel ('frequency(Hz)'); 
ylabel ('Power/Frequency(dB)');
grid on;
%%Modified Periodogram with Hamming Window. 
% Obtain the modified periodogram of an input EEG signal with no noise. The
% signal is 30001 samples in length. Obtain the modified periodogram using
% a Hamming window and default DFT length. The DFT length is the next power
% of two greater than the signal length, or 32786 points. 
%Because the signal is real-valued and has even length, the periodogram is one-sided and there are 32768/2+1 points.   
hamming_1=periodogram(eeg1,hamming(length(eeg1)));
subplot(2,2,3), 
plot (f,hamming_1);
title ('Periodogram using the hamming window');





z=inSdfsMat(10,:);
fs = 1000;

subplot(3,2,1)
plot(z)

L=length(z);
NEFT=2^nextpow2(L);
z_fft=abs(fft(z,NEFT));
freq=fs/2*linspace(0,1.5,NEFT/2+1);
subplot(3,2,4);
plot(freq,z_fft(1:length(freq)));

subplot(3,2,2)
z_ifft = abs(ifft([z_fft(1:10) zeros(1,1001-9)]));
plot(z_ifft,'r')
hold on
plot(diff(z_ifft,2))
%%
for ii = 1:size(inSdfsMat,1)
z=inSdfsMat(ii,stWin:enWin);
z_fft = abs(fft(z,length(z)));
z_ifft = abs(ifft(z_fft(1:50),length(z)));
plot(z) 
hold on 
plot(z_ifft)
st(ii) = find((z - z_ifft)>0,1);
line([st(ii) st(ii)],get(gca,'YLim'),'color','k')
hold off
pause
end 


%%


o=5;
wn=[1 10]*2/fs;
[b,a]=butter(o,wn,'bandpass');
figure;freqz(b,a,1000,fs);
figure(1)
[h,w]=freqz(b,a,1000,fs);
subplot(3,2,3);
plot(w,20*log10(abs(h)));

z_filt=filter(b,a,z);
subplot(3,2,5);
plot(z_filt,'b','linewidth', 2);


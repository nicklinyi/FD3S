

fprintf(1,'\nwavelet parameters:\n');
wavelet=input('wavelet: Gauss=1, Ricker=2, Dirac=3: ');
t_0=input('t_0: ');
d_t=input('d_t: ');
n_t=input('n_t: ');
t_s=input('t_s: ');
t_p=input('t_p: ');

% compute wavelet

if (wavelet==1)
    f=gauss(t_0,n_t,d_t,t_s,t_p);
elseif (wavelet==2)
    f=ricker(t_0,n_t,d_t,t_s,t_p);
elseif (wavelet==3)
    f=(1:(n_t+1))*0;
    f(2)=1;
end

t=t_0:d_t:n_t*d_t;

% plot signal

figure
subplot(2,1,1)
plot(t,f,'k');
xlabel('t in s');
title('wavelet, time domain');

k=[0:length(t)-1]/max(t);
ff=abs(fft(f));
subplot(2,1,2)
plot(k(1:round(end/2)),ff(1:round(end/2)),'k');
xlabel('\nu in Hz');
title('wavelet, frequency domain');


fid=fopen('../MODELS/stf','w');

for k=1:length(t)
    fprintf(fid,'%g\n',f(k));
end

fclose(fid);

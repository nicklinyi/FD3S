function f=ricker(t_0,n_t,d_t,t_s,t_p)

% function f=ricker(t_0,n_t,d_t,t_s,t_p)

%t_0=0;
%n_t=100;
%d_t=0.1;

%t_s=1.0;
%t_p=0.5;


t=t_0:d_t:n_t*d_t;

a=3.1415*(t-t_s).^2/t_p^2;
f=0.5*3.1415*(a-0.5).*exp(-a);

% plot
if 0
    figure
    subplot(2,1,1)
    plot(t,f,'k');
    xlabel('t in s');
    title('Ricker wavelet, time domain');

    k=[0:length(t)-1]/max(t);
    ff=abs(fft(f));
    subplot(2,1,2)
    plot(k(1:round(end/2)),ff(1:round(end/2)),'k');
    xlabel('\nu in Hz');
    title('Ricker wavelet, frequency domain');
end

function [f_center,t_max,F,T,tt]=z_dispersion(p,f_min,f_max,nf,trace)


% compute the Fourier transform of the spectrum

fp=z_fourier(p);
ffp=fp;

% considered frequency range

f_center=f_min:((f_max-f_min)/(nf-1)):f_max;
t=p.t0:p.dt:(p.nt-1)*p.dt+p.t0;
tt=zeros(nf,p.nt);
t_max=1:nf;

% create the Gaussian filter sequence

f=0:fp.df:(fp.nt-1)*fp.df;
index=1:nf;

for n=1:nf

    alpha=1700*f_center(n)^2;
    window=exp(-alpha*(f-f_center(n)).^2/f_center(n)^2);  
    ffp.spectra(trace,:)=fp.spectra(trace,:).*window;  
    iffp=z_ifourier(ffp); 
    
    tt(n,:)=abs(iffp.seismograms(trace,:));
    
    t_max(n)=find(max(abs(iffp.seismograms(trace,:)))==abs(iffp.seismograms(trace,:)));
    t_max(n)=t(t_max(n));
    
end

[F,T]=meshgrid(f_center,t);


figure
pcolor(F,T,tt')
shading interp
colormap bone

hold on
plot(f_center,t_max,'k')
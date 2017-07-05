function fp=z_fourier(p)

% function fp=z_fourier(p)
%
% computes the FOURIER spectrum of the seismograms contained in the
% seismogram structure p
%
% OUTPUT is a spectrum structure fp with the following components:
%
% fp.nr :       number of seismograms=number of spectra
% fp.nt :       number of points in both time and frequency
% fp.dt :       time increment of seismograms
% fp.t0 :       initial time of the seismograms
% fp.r_theta :  theta coordinate of the receivers
% fp.r_phi :    phi coordinate of the receivers
% fp.spectra :  spectra
% fp.df :       frequency increment


fp.nr=p.nr;
fp.nt=p.nt;
fp.dt=p.dt;
fp.t0=p.t0;
fp.r_theta=p.r_theta;
fp.r_phi=p.r_phi;

fp.spectra=zeros(p.nr, p.nt);


i=sqrt(-1);
KK=0:(p.nt-1);

for k=1:p.nr
    fp.spectra(k,:)=fft(p.seismograms(k,:));
    fp.spectra(k,:)=p.dt*fp.spectra(k,:).*exp((-2*pi*i*KK*p.t0)/(p.nt*p.dt));
end

% wrap around the high frequencies to negative frequencies
if 0
n=ceil(p.nt/2)
dummy=zeros([1 n]);

for k=1:p.nr
    dummy(1:n)=fp.spectra(k,1:n);
    for m=1:(n-1)
        fp.spectra(k,m)=fp.spectra(k,n+m)*exp(-2*pi*i*p.t0*p.dt);
    end
    fp.spectra(k,n:p.nt)=dummy(:);
end
end
fp.df=1/(p.nt*p.dt);
function ifp=z_ifourier(fp)

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


ifp.nr=fp.nr;
ifp.nt=fp.nt;
ifp.dt=fp.dt;
ifp.t0=fp.t0;
ifp.r_theta=fp.r_theta;
ifp.r_phi=fp.r_phi;

ifp.seismograms=zeros(fp.nr, fp.nt);
%size(fp.spectra)

i=sqrt(-1);
KK=0:(fp.nt-1);
%size(K)
for k=1:fp.nr
    fp.spectra(k,:)=fp.spectra(k,:).*exp((2*pi*i*KK*fp.t0)/(fp.nt*fp.dt));
    ifp.seismograms(k,:)=ifft(fp.spectra(k,:));
    ifp.seismograms(k,:)=ifp.seismograms(k,:)/fp.dt;
end

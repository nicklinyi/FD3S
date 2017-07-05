function cp=z_cut_seismograms(p,min_time,max_time)

% function cp=z_cut_seismograms(p,min_time,max_time)
%
% Cuts seismograms on the time axis and assigns new initial times to the
% seismogram structure
%
% INPUT:    
%
% p:        seismogram structure
% min_time: minimum time
% max_time: maximum time
%
% OUTPUT:
% 
% cp:       modified seismogram structure
%
% last modified: 30 June, 2006

cp.nr=p.nr;
cp.dt=p.dt;
cp.r_theta=p.r_theta;
cp.r_phi=p.r_phi;
cp.t0=min_time;

t=p.t0:p.dt:(p.nt-1)*p.dt+p.t0;

% determine startin and ending indeces

nmin=round((min_time-p.t0)/p.dt);
nmax=round((max_time-p.t0)/p.dt);

cp.nt=nmax-nmin+1;
cp.seimograms=zeros(p.nr,cp.nt);

for k=1:p.nr
    cp.seismograms(k,:)=p.seismograms(k,nmin:nmax);
end
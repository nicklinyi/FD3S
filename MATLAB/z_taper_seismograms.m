function tp=z_taper_seismograms(p,t_min_array,t_max_array,width)

% function tp=taper_seismograms
% (
%   p,              | seismogram structure
%   t_min_array,    | minimum taper boundaries
%   t_max_array,    | maximum taper boundaries
%   width           | width of the cosine part
% )
%
% tapers seismograms
% output is a structure containing the tapered seismograms


tp.nr=p.nr;
tp.nt=p.nt;
tp.dt=p.dt;
tp.r_theta=p.r_theta;
tp.r_phi=p.r_phi;
tp.t0=p.t0;

window=ones(p.nr,p.nt);
t=p.t0:p.dt:p.dt*(p.nt-1)+p.t0;

for k=1:p.nr
    window(k,:)=window(k,:).*(t>t_min_array(k)).*(t<t_max_array(k));
    window(k,:)=window(k,:)+(0.5+0.5*cos(pi*(t_max_array(k)-t)/(width))).*(t>=t_max_array(k)).*(t<t_max_array(k)+width);
    window(k,:)=window(k,:)+(0.5+0.5*cos(pi*(t_min_array(k)-t)/(width))).*(t>t_min_array(k)-width).*(t<=t_min_array(k));
    tp.seismograms(k,:)=p.seismograms(k,:).*window(k,:);
end
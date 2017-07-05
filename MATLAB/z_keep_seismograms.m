function kp=z_keep_seismograms(p,keep_array)

kp.nr=length(keep_array);
kp.nt=p.nt;
kp.dt=p.dt;
kp.r_theta=p.r_theta(keep_array);
kp.r_phi=p.r_phi(keep_array);
kp.seismograms=p.seismograms(keep_array,:);
kp.t0=p.t0;
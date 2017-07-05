function z_plot_seismograms(p,array,spacing)

% function z_plot_seismograms(p,array,spacing)
%
% plots seismograms from a seismogram structure
%
% INPUT:    
%
% p:        seismogram structure
% array:    array of the indeces of the seismograms that will be plotted,
%           enter [] if all seismograms are to be plotted
% spacing:  spacing of seismograms in the plot as multiple of the maximum
%           amplitude. A reasonable value is 2.0 .
%
% OUTPUT:   none
%
% last modified: 30 June, 2006

m=max(max(abs(p.seismograms)));
t=p.t0:p.dt:((p.nt-1)*p.dt)+p.t0;

figure
hold on

if (length(array)==0)

    for ir=1:p.nr
        plot(t,(ir-1)*spacing*m+p.seismograms(ir,:),'b');
        text(1.05*t(end),(ir-1)*spacing*m,num2str(ir));
        text(-0.09*t(end)+p.t0,(ir-1)*spacing*m+0.1*m,num2str(p.r_theta(ir)));
        text(-0.01*t(end)+p.t0,(ir-1)*spacing*m+0.1*m,num2str(p.r_phi(ir)));
    end

    axis([-0.1*t(end)+p.t0 1.1*t(end) -1.1*m (p.nr-1)*spacing*m+1.1*m]);
    grid on
    
else
    
    for ir=1:length(array)
        plot(t,(ir-1)*spacing*m+p.seismograms(array(ir),:),'b');
        text(1.05*t(end)+5,(ir-1)*spacing*m,num2str(ir));
        text(-0.09*t(end)+p.t0,(ir-1)*spacing*m+0.1*m,num2str(p.r_theta(array(ir))));
        text(-0.01*t(end)+p.t0,(ir-1)*spacing*m+0.1*m,num2str(p.r_phi(array(ir))));
    end

    axis([-0.1*t(end)+p.t0 1.1*t(end) -1.1*m (length(array)-1)*spacing*m+1.1*m]);
    grid on
    
end

xlabel('time [s]');
ylabel('velocity [m/s]');

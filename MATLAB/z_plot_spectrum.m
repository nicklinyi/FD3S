function z_plot_spectrum(fp,array,spacing)

% function z_plot_spectrum(fp,array,spacing)
%
% plots the left part of the spectrum computed with
% z_fourier

m=max(max(abs(fp.spectra)));

f=0:fp.df:(fp.nt-1)*fp.df;

n=ceil(fp.nt/2);
%f=-(0.5*fp.nt*fp.df):fp.df:(0.5*fp.nt*fp.df)-fp.df;

figure
hold on

x=abs(fp.spectra);

if (length(array)==0)

    for ir=1:fp.nr
        plot(f,(ir-1)*spacing*m+x(ir,:),'b');
        text(1.1*f(end),(ir-1)*spacing*m+0.1*m,num2str(ir));
        text(-0.2*f(end),(ir-1)*spacing*m+0.1*m,num2str(fp.r_theta(ir)));
        text(-0.1*f(end),(ir-1)*spacing*m+0.1*m,num2str(fp.r_phi(ir)));
    end

    axis([-0.2*f(end) 1.2*f(end) -0.1*m fp.nr*spacing*m+0.1*m]);
    grid on
    
else
    
    for ir=1:length(array)
        plot(f,(ir-1)*spacing*m+x(array(ir),:),'b');
        text(1.1*f(end),(ir-1)*spacing*m+0.1*m,num2str(ir));
        text(-0.2*f(end),(ir-1)*spacing*m+0.1*m,num2str(fp.r_theta(array(ir))));
        text(-0.1*f(end),(ir-1)*spacing*m+0.1*m,num2str(fp.r_phi(array(ir))));
    end

    axis([-0.2*f(end) 1.2*f(end) -0.1*m fp.nr*spacing*m+0.1*m]);
    grid on
    
end

xlabel('frequency [Hz]');
ylabel('amplitude spectrum');

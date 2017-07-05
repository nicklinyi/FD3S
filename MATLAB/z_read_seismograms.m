function p=z_read_seismograms(array,spacing)

% function p=z_read_seismograms(array,spacing)
%
% reads seismograms out of the output files and plots them
%
% INPUT:
%
% array:    array of the indeces of the seismograms that will be plotted,
%           enter [] if all seismograms are to be plotted
%
% spacing:  spacing of seismograms in the plot as multiple of the maximum
%           amplitude. A reasonable value is 2.0 .
%
% OUTPUT:   a seismogram structure
%
% last modified: 30 June, 2006

component=input('component (theta, phi, r): ','s');
processor=input('processor: ','s');

filename=['../DATA/OUTPUT/seismogram_' component processor];

fid=fopen(filename,'r');

fgets(fid);
fscanf(fid,'%c',4);
p.nr=fscanf(fid,'%d',1);
fscanf(fid,'%c',5);
p.nt=fscanf(fid,'%d',1);
fscanf(fid,'%c',5);
p.dt=fscanf(fid,'%g',1);

p.r_theta=zeros(1,p.nr);
p.r_phi=zeros(1,p.nr);
p.seismograms=zeros(p.nr,p.nt);

p.t0=0.0;

for ir=1:p.nr
    fscanf(fid,'%c',9);
    p.r_theta(ir)=fscanf(fid,'%g',1);
    p.r_theta(ir)=p.r_theta(ir)*180/pi;
    fscanf(fid,'%c',6);
    p.r_phi(ir)=fscanf(fid,'%g',1);
    p.r_phi(ir)=p.r_phi(ir)*180/pi;
    for it=1:p.nt    
        p.seismograms(ir,it)=fscanf(fid,'%g',1);
    end
end

% plot seismograms

m=max(max(abs(p.seismograms)));
t=0:p.dt:((p.nt-1)*p.dt);

figure
hold on

if (length(array)==0)

    for ir=1:p.nr
        plot(t,ir*spacing*m+p.seismograms(ir,:),'b');
        text(1.05*t(end),ir*spacing*m,num2str(ir));
        text(-0.09*t(end),ir*spacing*m+0.1*m,num2str(p.r_theta(ir)));
        text(-0.01*t(end),ir*spacing*m+0.1*m,num2str(p.r_phi(ir)));
    end

    axis([-0.1*t(end) 1.1*t(end) 0 p.nr*spacing*m+1.1*m]);
    grid on
    
else
    
    for ir=1:length(array)
        plot(t,ir*spacing*m+p.seismograms(array(ir),:),'b');
        text(1.05*t(end)+5,ir*spacing*m,num2str(ir));
        text(-0.09*t(end),ir*spacing*m+0.1*m,num2str(p.r_theta(array(ir))));
        text(-0.01*t(end),ir*spacing*m+0.1*m,num2str(p.r_phi(array(ir))));
    end

    axis([-0.1*t(end) 1.1*t(end) 0 p.nr*spacing*m+1.1*m]);
    grid on
    
end

xlabel('time [s]');
ylabel('velocity [m/s]');

fclose(fid);

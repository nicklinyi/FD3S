%==========================================================================
% read boxfile
%==========================================================================

fpt=fopen('../MODELS/boxfile','r');

for i=1:14
    fgets(fpt);
end

p.pp=fscanf(fpt,'%i',1);
p.px=fscanf(fpt,'%i',1);
p.py=fscanf(fpt,'%i',1);
p.pz=fscanf(fpt,'%i',1);

fgets(fpt);
fgets(fpt);

p.mi=[1:p.pp];
p.mi_x=[1:p.pp];
p.mi_y=[1:p.pp];
p.mi_z=[1:p.pp];

p.ix_min=[1:p.pp];
p.ix_max=[1:p.pp];
p.iy_min=[1:p.pp];
p.iy_max=[1:p.pp];
p.iz_min=[1:p.pp];
p.iz_max=[1:p.pp];

p.x_min=[1:p.pp];
p.x_max=[1:p.pp];
p.y_min=[1:p.pp];
p.y_max=[1:p.pp];
p.z_min=[1:p.pp];
p.z_max=[1:p.pp];

for i=1:p.pp
    
    p.mi(i)=fscanf(fpt,'%d',1);
    p.mi_x(i)=fscanf(fpt,'%d',1);
    p.mi_y(i)=fscanf(fpt,'%d',1);
    p.mi_z(i)=fscanf(fpt,'%d',1);
    p.ix_min(i)=fscanf(fpt,'%i',1);
    p.ix_max(i)=fscanf(fpt,'%i',1);
    p.iy_min(i)=fscanf(fpt,'%i',1);
    p.iy_max(i)=fscanf(fpt,'%i',1);
    p.iz_min(i)=fscanf(fpt,'%i',1);
    p.iz_max(i)=fscanf(fpt,'%i',1);
    
    p.x_min(i)=fscanf(fpt,'%lf',1);
    p.x_max(i)=fscanf(fpt,'%lf',1);
    p.y_min(i)=fscanf(fpt,'%lf',1);
    p.y_max(i)=fscanf(fpt,'%lf',1);
    p.z_min(i)=fscanf(fpt,'%lf',1);
    p.z_max(i)=fscanf(fpt,'%lf',1);
    
    fgets(fpt);
    fgets(fpt);
    
end

fclose(fpt);

%==========================================================================
% determine filename
%==========================================================================

n=input('processor: ');

filename_lambda=['../MODELS/lambda' num2str(n)];
filename_mu=['../MODELS/mu' num2str(n)];
filename_rhoinv=['../MODELS/rhoinv' num2str(n)];

fpt_lambda=fopen(filename_lambda,'r');
fpt_mu=fopen(filename_mu,'r');
fpt_rhoinv=fopen(filename_rhoinv,'r');

fprintf(1,'opening file %s\n',filename_lambda);
fprintf(1,'opening file %s\n',filename_mu);
fprintf(1,'opening file %s\n',filename_rhoinv);


%==========================================================================
% read profiles
%==========================================================================

nz=p.iz_max(n+1)-p.iz_min(n+1);
ny=p.iy_max(n+1)-p.iy_min(n+1);
nx=p.ix_max(n+1)-p.ix_min(n+1);

lambda=zeros([1 nz+1]);
mu=zeros([1 nz+1]);
rhoinv=zeros([1 nz+1]);

for k=1:nz+1
    for j=1:ny+1
        for i=1:nx+1
            lambda(k)=fscanf(fpt_lambda,'%lg',1);
            mu(k)=fscanf(fpt_mu,'%lg',1);
            rhoinv(k)=fscanf(fpt_rhoinv,'%lg',1);
        end
    end
    fprintf(1,'%d ',k);
end

%==========================================================================
% plot output 
%==========================================================================

dz=(p.z_max(n+1)-p.z_min(n+1))/nz;
z=p.z_max(n+1):-dz:p.z_min(n+1);

figure
plot(z,mu,'k');
hold on
plot(z,lambda,'r');
title('\mu (black), \lambda (red)');

cp=sqrt((2*mu+lambda).*rhoinv);
cs=sqrt(mu.*rhoinv);
figure
plot(z,cp,'k');
hold on
plot(z,cs,'r');
title('c_p (black), c_s (red)');



%==========================================================================
% clean up
%==========================================================================

fclose(fpt_lambda);
fclose(fpt_mu);
fclose(fpt_rhoinv);
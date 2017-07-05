%==========================================================================
% read boxfile
%==========================================================================

figure

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
plane=input('plane (theta, phi, r): ','s');
comp=input('component (theta, phi, r): ', 's');
it=input('iteration: ');

filename=['../DATA/OUTPUT/vad_' comp '_' plane 'plane' num2str(n) '_' num2str(it)];
fpt=fopen(filename,'r');

fprintf(1,'opening file %s\n',filename);

if (strcmp(plane,'theta'))
    
    field=zeros(p.iy_max(n+1)-p.iy_min(n+1)+1,p.iz_max(n+1)-p.iz_min(n+1)+1);
    
    for k=0:1:p.iz_max(n+1)-p.iz_min(n+1)
        for j=0:1:p.iy_max(n+1)-p.iy_min(n+1)
            field(j+1,k+1)=fscanf(fpt,'%lf',1);
            %fprintf(1,'.',1);
        end
    end
    
    y=(p.iy_min(n+1)-2):(p.iy_max(n+1)+2);
    z=(p.iz_min(n+1)-2):(p.iz_max(n+1)+2);
    %image(y,z,field');
    %caxis([-max(max(abs(field))) max(max(abs(field)))]);
    pcolor(field')
    
end

if (strcmp(plane,'phi'))
    
    field=zeros(p.ix_max(n+1)-p.ix_min(n+1)+1,p.iz_max(n+1)-p.iz_min(n+1)+1);
    
    for k=0:1:p.iz_max(n+1)-p.iz_min(n+1)
        for i=0:1:p.ix_max(n+1)-p.ix_min(n+1)
            field(i+1,k+1)=fscanf(fpt,'%lf',1);
            %fprintf(1,'.',1);
        end
    end
    
    %x=(p.ix_min(n+1)-2):(p.ix_max(n+1)+2);
    %z=(p.iz_min(n+1)-2):(p.iz_max(n+1)+2);
    %image(x,z,field');
    pcolor(field')
    
end

if (strcmp(plane,'r'))
    
    field=zeros(p.ix_max(n+1)-p.ix_min(n+1)+1,p.iy_max(n+1)-p.iy_min(n+1)+1);
    
    for j=0:1:p.iy_max(n+1)-p.iy_min(n+1)
        for i=0:1:p.ix_max(n+1)-p.ix_min(n+1)
            field(i+1,j+1)=fscanf(fpt,'%lf',1);
            %fprintf(1,'.',1);
        end
    end
    
    %x=(p.ix_min(n+1)-2):(p.ix_max(n+1)+2);
    %y=(p.iy_min(n+1)-2):(p.iy_max(n+1)+2);
    %image(x,y,field');
    pcolor(field')
    
end


fclose(fpt);
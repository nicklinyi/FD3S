%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots radial profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read model parameters

fid=fopen('model_parameters','r');

nx=fscanf(fid,'%i',1); fgets(fid);
ny=fscanf(fid,'%i',1); fgets(fid);
nz=fscanf(fid,'%i',1); fgets(fid);
z_max=fscanf(fid,'%lf',1); fgets(fid);
z_min=fscanf(fid,'%lf',1); fgets(fid);

fclose(fid);

% read radial profile

profile=zeros(1,nz+1);
z=zeros(1,nz+1);

filename=input('filename: ','s');
fid=fopen(filename,'r');

for i=0:nz
    for j=0:ny
        for k=0:nx
            profile(i+1)=fscanf(fid,'%lf',1);
        end
    end
end

fclose(fid);

% plot results

dz=(z_max-z_min)/nz;
for k=0:nz
    z(k+1)=z_max-dz*k;
end

figure
plot(z,profile,'k');
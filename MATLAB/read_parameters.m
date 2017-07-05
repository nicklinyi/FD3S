% copy Par-file to MATLAB directory
cd ..
!cp Par MATLAB
!cp params.h MATLAB
cd MATLAB


% read file

par_file=fopen('Par','r');

dummy=fscanf(par_file,'%c',79);
p.seisfile=fscanf(par_file,'%s',1); fgets(par_file); 
p.recfile=fscanf(par_file,'%s',1); fgets(par_file); 
p.nt=fscanf(par_file,'%i',1); fgets(par_file); 
p.dt=fscanf(par_file,'%lf',1); fgets(par_file); 
p.izfree=fscanf(par_file,'%i',1); fgets(par_file);
fgets(par_file);
p.ibound=fscanf(par_file,'%i',1); fgets(par_file); 
p.nabs=fscanf(par_file,'%i',1); fgets(par_file); 
fgets(par_file);
p.xxs=fscanf(par_file,'%lf',1); fgets(par_file); 
p.yys=fscanf(par_file,'%lf',1); fgets(par_file); 
p.zzs=fscanf(par_file,'%lf',1); fgets(par_file); 
p.srctype=fscanf(par_file,'%i',1); fgets(par_file); 
p.dip=fscanf(par_file,'%lf',1); fgets(par_file); 
p.rake=fscanf(par_file,'%lf',1); fgets(par_file); 
p.strike=fscanf(par_file,'%lf',1); fgets(par_file); 
fgets(par_file);
%p.isamp=fscanf(par_file,'%i',1); fgets(par_file); 
p.ssamp=fscanf(par_file,'%i',1); fgets(par_file); 
p.iflush=fscanf(par_file,'%i',1); fgets(par_file); 
p.iplanex=fscanf(par_file,'%i',1); fgets(par_file); 
p.iplaney=fscanf(par_file,'%i',1); fgets(par_file); 
p.iplanez=fscanf(par_file,'%i',1); fgets(par_file); 
p.icheck=fscanf(par_file,'%i',1); fgets(par_file); 
p.xmin=fscanf(par_file,'%lf',1); fgets(par_file); 
p.xmax=fscanf(par_file,'%lf',1); fgets(par_file); 
p.ymin=fscanf(par_file,'%lf',1); fgets(par_file); 
p.ymax=fscanf(par_file,'%lf',1); fgets(par_file); 
p.zmin=fscanf(par_file,'%lf',1); fgets(par_file); 
p.zmax=fscanf(par_file,'%lf',1); fgets(par_file); 
p.rotdiv=fscanf(par_file,'%i',1); fgets(par_file); 
p.fs_model=fscanf(par_file,'%i',1); fgets(par_file); 
p.elast_inter=fscanf(par_file,'%i',1); fgets(par_file);
p.is_aniso=fscanf(par_file,'%i',1); fgets(par_file);
p.is_diss=fscanf(par_file,'%i',1);

fclose(par_file);
par_file=fopen('params.h','r');

read_dummy(par_file);
p.nx=fscanf(par_file,'%i',1); read_dummy(par_file);
p.ny=fscanf(par_file,'%i',1); read_dummy(par_file);
p.nz=fscanf(par_file,'%i',1); read_dummy(par_file);
p.maxnt=fscanf(par_file,'%i',1); read_dummy(par_file);
p.maxnr=fscanf(par_file,'%i',1); read_dummy(par_file);

fclose(par_file);


% clean up

!rm params.h
!rm Par
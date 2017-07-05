function convp=z_convolve(p)

% z_convolve assumes that the function with which the seismograms are
% convolved, is causal, i.e., zero prior to t=0. Hence, one should choose
% the parameters for the function reasonably.

convp=p;
t=p.t0:p.dt:(p.nt-1)*p.dt+p.t0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelet generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavelet=input('wavelet: (1=Gauss, 2=Ricker)');

if (wavelet==1)
    t_s=input('t_s= ');
    t_p=input('t_p= ');
elseif (wavelet==2)
    t_s=input('t_s= ');
    t_p=input('t_p= ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (wavelet==1)
    
    for k=1:p.nr
        for n=1:p.nt
            convp.seismograms(k,n)=sum(p.seismograms(k,:).*fgauss(p.t0+(n-1)*p.dt-t,t_s,t_p));
        end
        fprintf(1,'.');
    end
    fprintf(1,'\n');

end

if (wavelet==2)
    
    for k=1:p.nr
        for n=1:p.nt
            convp.seismograms(k,n)=sum(p.seismograms(k,:).*fricker(p.t0+(n-1)*p.dt-t,t_s,t_p));
        end
        fprintf(1,'.');
    end
    fprintf(1,'\n');

end
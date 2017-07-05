function f=fgauss(t,t_s,t_p)

f=exp(-(t-t_p).^2/t_s^2);
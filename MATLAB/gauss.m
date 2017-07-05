function f=gauss(t_0,n_t,d_t,t_s,t_p)

t=t_0:d_t:n_t*d_t;
f=exp(-(t-t_p).^2/t_s^2);
function f=fricker(t,t_s,t_p)

f=-2*(t-t_p).*exp(-(t-t_p).^2/t_s^2)/t_s^2;
function y_n = rk4(fn,y,u,ts,i)

t0=0;
t=t0+(i-1)*ts;
h=ts;

k1=h*feval(fn,y,t,u);
k2=h*feval(fn,y+k1/2,t+h/2,u);
k3=h*feval(fn,y+k2/2,t+h/2,u);
k4=h*feval(fn,y+k3,t+h,u);

y_n=y+(k1+2*k2+2*k3+k4)/6;
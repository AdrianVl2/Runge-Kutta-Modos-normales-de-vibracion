k1=10; k2=10; k3=10; m1=1; m2=1;
%F1 = @(t) 5*sin(2*t);   % fuerza en m1
%F2 = @(t) 0*t;          % sin fuerza en m2
[t,x1,v1,x2,v2] = rk4_segundo_orden_masas(k1,k2,k3,m1,m2,0,20,0.0001,...
                                          1,0, 0,0);
plot(t,x1,t,x2); legend('x_1','x_2'), xlabel t, ylabel desplazamiento

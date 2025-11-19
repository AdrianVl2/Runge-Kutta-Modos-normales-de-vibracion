function [t,x1,v1,x2,v2] = rk4_segundo_orden_masas(k1,k2,k3,m1,m2,t0,tf,h,...
                                                   x10,v10,x20,v20)
% RK4 directo para sistema masa-resorte acoplado de 2º orden con forzamiento.
% F1,F2: function handles @(t) ...
% Condiciones iniciales: x10,v10,x20,v20

t = t0:h:tf;  N = numel(t);
x1 = zeros(1,N); v1 = zeros(1,N);
x2 = zeros(1,N); v2 = zeros(1,N);

x1(1)=x10; v1(1)=v10;  x2(1)=x20; v2(1)=v20;

a1 = @(tt,X1,V1,X2,V2) (-(k1+k2)*X1 + k2*X2 )/m1;
a2 = @(tt,X1,V1,X2,V2) ( k2*X1 - (k2+k3)*X2)/m2;

for n=1:N-1
    tn=t(n);

    % k1
    k11 = h*a1(tn,x1(n),v1(n),x2(n),v2(n));
    k21 = h*a2(tn,x1(n),v1(n),x2(n),v2(n));

    % mid 1 (con k11,k21)
    x1m = x1(n) + (h/2)*v1(n) + (h^2/8)*k11;
    v1m = v1(n) + k11/2;
    x2m = x2(n) + (h/2)*v2(n) + (h^2/8)*k21;
    v2m = v2(n) + k21/2;

    % k2
    k12 = h*a1(tn+h/2,x1m,v1m,x2m,v2m);
    k22 = h*a2(tn+h/2,x1m,v1m,x2m,v2m);

    % mid 2 (con k12,k22)
    x1m = x1(n) + (h/2)*(v1(n)+k11/2) + (h^2/8)*k12;
    v1m = v1(n) + k12/2;
    x2m = x2(n) + (h/2)*(v2(n)+k21/2) + (h^2/8)*k22;
    v2m = v2(n) + k22/2;

    % k3
    k13 = h*a1(tn+h/2,x1m,v1m,x2m,v2m);
    k23 = h*a2(tn+h/2,x1m,v1m,x2m,v2m);

    % end (con k13,k23)
    x1e = x1(n) + h*(v1(n)+k12/2) + (h^2/2)*k13;
    v1e = v1(n) + k13;
    x2e = x2(n) + h*(v2(n)+k22/2) + (h^2/2)*k23;
    v2e = v2(n) + k23;

    % k4
    k14 = h*a1(tn+h,x1e,v1e,x2e,v2e);
    k24 = h*a2(tn+h,x1e,v1e,x2e,v2e);

    % actualizar
    v1(n+1) = v1(n) + (k11 + 2*k12 + 2*k13 + k14)/6;
    v2(n+1) = v2(n) + (k21 + 2*k22 + 2*k23 + k24)/6;

    x1(n+1) = x1(n) + h*( v1(n) + (h/6)*(k11+k12+k13) );
    x2(n+1) = x2(n) + h*( v2(n) + (h/6)*(k21+k22+k23) );
end
end

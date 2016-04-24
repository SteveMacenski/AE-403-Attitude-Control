function hw5_prob1

clc;
%% theta(t) v(t) and u(t) 1d
wn = [10^-2 10^-1 10^0 10^1 10^2];
J = 10;
Jw = 1; 
A = [0 1; 0 0];
B = [0 1/J]';
x0 = [0.1 0]';
v0 = 0;
for j=1:length(wn)

    t = linspace(0,9,400);

    k = [wn(j)^2*J 2*J/sqrt(2)*wn(j)];
    
    for i=1:length(t)
        x(:,i) = expm((A-B*k)*t(i))*x0;
    end
    figure(1) %theta(t)
    plot(t,x(1,:))
    hold on;
    
    figure(2) %u(t)
    plot(t,-k*x)
    hold on;

    for i=1:length(t)
        v(:,i) = -k*(1/Jw)*inv(A-B*k)*(expm((A-B*k)*t(i))-eye(2))*x0;
    end
    
    figure(3) %v(t)
    plot(t,v)
    hold on;
end
figure(1)
legend('theta wn=10e-2','theta wn=10e-1','theta wn=10e0','theta wn=10e1','theta wn=10e2')
title('\theta(t) - Steve Macenski - P1.d')

figure(2)
legend('u(t) wn=10e-2','u(t) wn=10e-1','u(t) wn=10e0','u(t) wn=10e1','u(t) wn=10e2')
title('u(t) - Steve Macenski - P1.d')

figure(3)
legend('v wn=10e-2','v wn=10e-1','v wn=10e0','v wn=10e1','v wn=10e2')
title('v(t) - Steve Macenski - P1.d')

%% theta(t) v(t) u(t) for 1e: using LQR 
R = [10^-2 10^-1 10^0 10^1 10^2];
J = 10;
Jw = 1;
A = [0 1; 0 0];
B = [0 1/J]'; 
x0 = [0.1 0]';
v0 = 0;
Q = eye(2);
for j=1:length(R)

    t = linspace(0,9,1400);

    k = lqr(A,B,Q,R(j));

    for i=1:length(t)
        x(:,i) = expm((A-B*k)*t(i))*x0;
    end
    figure(4) %theta(t)
    plot(t,x(1,:))
    hold on;
    
    figure(5) %u(t)
    plot(t,-k*x)
    hold on;
    
    for i=1:length(t)
        v(:,i) = -k*(1/Jw)*inv(A-B*k)*(expm((A-B*k)*t(i))-eye(2))*x0;
    end
    
    figure(6) %v(t)
    plot(t,v)
    hold on;
end 
figure(4)
legend('theta r=10e-2','theta r=10e-1','theta r=10e0','theta r=10e1','theta r=10e2')
title('\theta(t) - Steve Macenski - P1.e')

figure(5)
legend('u(t) r=10e-2','u(t) r=10e-1','u(t) r=10e0','u(t) r=10e1','u(t) r=10e2')
title('u(t) - Steve Macenski - P1.e')

figure(6)
legend('v(t) r=10e-2','v(t) r=10e-1','v(t) r=10e0','v(t) r=10e1','v(t) r=10e2')
title('v(t) - Steve Macenski - P1.e')



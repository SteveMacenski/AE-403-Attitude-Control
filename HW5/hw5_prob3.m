function hw5_prob3

clc;
%% x(t) problem 3
J1 = 12;
J2 = 14;
J3 = 8; 
Jw = 1; 
n = .0011;


A = [      0                 0               1               0          0     0;...
           0                 0               0               1          0     0;...
     -4*n^2*(J2-J3)/J1       0               0         n*(1-(J2-J3)/J1) 0     0;...
           0           n^2*(J1-J2)/J3 -n*(1+(J1-J2)/J3)      0          0     0;...
           0                 0               0               0          0     n;...
           0                 0               0               0         -n     0];
       
B = [0     0;...
     0     0;...
     1/J1  0;...
     0   1/J3;...
     -1/Jw 0;...
     0  -1/Jw];
 
rank(ctrb(A,B)); %6 is full rank, this full rank 
x0 = [0.1 0.5 0 0 0 0]';

Q = eye(6);
R = 10^0.*eye(2); 

t = linspace(0,300,300);

k = lqr(A,B,Q,R)


for i=1:length(t)
    x(:,i) = expm((A-B*k)*t(i))*x0;
end

% figure(1) %x(t)
% plot(t,x)
% hold on;
% legend('x1(t) - \theta_1','x2(t) - \theta_3','x3(t) - \theta_1dot','x4(t) - \theta_3dot','x5(t) - v_1','x6(t) - v_3')
% title('x(t) - Steve Macenski - P3.c');

figure(2) 
plot(t,x(1,:));
title('x1(t) - Steve Macenski - P3.c');
legend('x1(t) - \theta_1');

figure(3) 
plot(t,x(2,:));
title('x2(t) - Steve Macenski - P3.c');
legend('x2(t) - \theta_3');

figure(4) 
plot(t,x(3,:));
title('x3(t) - Steve Macenski - P3.c');
legend('x3(t) - \theta_1dot');

figure(5) 
plot(t,x(4,:));
title('x4(t) - Steve Macenski - P3.c');
legend('x4(t) - \theta_3dot');

figure(6) 
plot(t,x(5,:));
title('x5(t) - Steve Macenski - P3.c');
legend('x5(t) - v_1');

figure(7) 
plot(t,x(6,:));
title('x6(t) - Steve Macenski - P3.c');
legend('x6(t) - v_3');




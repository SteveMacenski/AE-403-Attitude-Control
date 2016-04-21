function hw3_prob5
clf;
woz = .0010471975; %rad/s
n = 10.4719; %rad/s
x0 = [0;0;n;woz]; %initial spin conditions [w1 w2 w3 v]
J = [4.83 0 0; 0 4.83 0; 0 0 10]; %moment of inertia 
Jd = 1; 
c = 1;
tMax = 500; 
tRate = 1;
t = linspace(0,tMax,tMax*tRate);

[t,x] = ode45(@(t,x) f(t,x,J,Jd,c),t,x0);  

figure(1) 
plot(t,x(:,1)); 
title('Steve Macenski')
ylabel('w1 (rad/s)')
xlabel('time (s)')

figure(2) 
plot(t,x(:,2)); 
title('Steve Macenski')
ylabel('w2 (rad/s)')
xlabel('time (s)')

figure(3) 
plot(t,x(:,3)); 
title('Steve Macenski')
ylabel('w3 (rad/s)')
xlabel('time (s)')

figure(4) 
plot(t,x(:,4)); 
title('Steve Macenski')
ylabel('v (rad/s)')
xlabel('time (s)')


function xdot = f(t,x,J,Jd,c)
w = x; %     [w1 w2 w3 v]
Jt = J(1,1);
Ja = J(3,3);
v = w(4);
wdot = [((Ja-Jt)*w(2)*w(3)-Jd*w(3)*v)/-Jt; (1/(Jt-Jd))*(w(1)*w(3)*(Ja-Jt)+c*v); (-Jd*w(1)*v)/Ja; (-c*v/Jd)-((1/(Jt-Jd))*(w(1)*w(3)*(Ja-Jt)+c*v))];

%fprintf('time: %f\n',t);
 xdot = wdot;
 


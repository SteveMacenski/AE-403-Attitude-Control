function hw5_prob2
% animation video: https://youtu.be/GnYS6B8IOGI
clc;
%% theta(t) v(t) and u(t) 2
J1 = 12;
J2 = 14;
J3 = 8; 
Jw = 1; 
n = .0011;

A = [0 1 0; 3*n^2*(J3-J1)/J2 0 0; 0 0 0];
B = [0 1/J2 -1/Jw]';
x0 = [0.5 0 0]';

Q = eye(3);
R = 10^0; 

t = linspace(0,160,300);

k = lqr(A,B,Q,R);

for i=1:length(t)
    x(:,i) = expm((A-B*k)*t(i))*x0;
end

figure(1) %x(t)
plot(t,x)
hold on;


figure(1)
legend('x1(t) - \theta_2','x2(t) - \theta_2dot','x3(t) - v_2')
title('x(t) - Steve Macenski - P2.c')




%% Animation and integration stuff 

%my movie stuff 
params.movie_filename = 'Steve_2_animation.avi'; % names

myV = VideoWriter(params.movie_filename);
myV.Quality = 100; %quality is the quality of video taken - keep at 100, 80 or less is just garbage
open(myV); %opens the video for inputting frames


% INITIAL CONDITIONS
R0 = [1 0 0;0 1 0; 0 0 1]; %initially aligned so it is an identity matrix
J = [J1 0 0; 0 J2 0; 0 0 J3];
tMax = 90; 
% DEFINE THE TIME INTERVAL
tRate = 30;
t = linspace(0,tMax,tMax*tRate);

% INTEGRATE THE ANGULAR RATE EQUATIONS **AND** EULER'S EQUATIONS
[t,x] = ode45(@(t,x) f(t,x,J1,J2,J3,Jw,n,A,B,x0,Q,R,k),t,[RtoX(R0); x0]);

% CREATE A BOX AND SOME AXES
Hb = 1;
Mb = 6*(J(1,1)+J(2,2)-J(3,3))/(Hb^2);
Wb = sqrt((6/Mb)*(J(1,1)-J(2,2)+J(3,3)));
Lb = sqrt((6/Mb)*(-J(1,1)+J(2,2)+J(3,3)));
[p1,faces,colors] = makebox(Lb,Wb,Hb);
pAxis = [0 0; 0 0; 0 1];
pAngVel = [zeros(3,1) x0/norm(x0)];
pAngMom = [zeros(3,1) J*x0/norm(J*x0)];

% SETUP THE PLOT
figure(2);
clf;
axis(1.25*[-1 1 -1 1 -1 1]);
axis equal;
hold on;
plot3(0,0,0,'k.','markersize',16);
hAxis = line(pAxis(1,:),pAxis(2,:),pAxis(3,:));
set(hAxis,'linewidth',2,'color','r');
hAngVel = line(pAngVel(1,:),pAngVel(2,:),pAngVel(3,:));
set(hAngVel,'linewidth',2,'color','g');
hAngMom = line(pAngMom(1,:),pAngMom(2,:),pAngMom(3,:));
set(hAngMom,'linewidth',2,'color','b');
hBox = patch('Vertices',p1','Faces',faces,...
          'CData',colors,'FaceColor','flat');
hTitle = title(sprintf('t = %4.2f',0));
lighting flat
light('Position',[0 -2 -1])
light('Position',[0 -2 1])
xlabel('x');
ylabel('y');
zlabel('z');

% ANIMATE THE RESULTS
w = x(:,10:end);
x = x(:,1:9);
firsttime = 1;
i = 1;
dt = max(t)/(length(t)-1);
tic;
while (i<length(t))
    if (toc > dt)
        tic;
        i = i+1;
        R = XtoR(x(i,:));
        p0 = R*p1;
        set(hBox,'Vertices',p0');
        pAxis0 = R*pAxis;
        pAngVel = [zeros(3,1) w(i,:)'/norm(w(i,:)')];
        pAngVel0 = R*pAngVel;
        pAngMom = [zeros(3,1) J*(w(i,:)')/norm(J*(w(i,:)'))];
        pAngMom0 = R*pAngMom;
        set(hAxis,'xdata',pAxis0(1,:),'ydata',pAxis0(2,:),'zdata',pAxis0(3,:));
        set(hAngVel,'xdata',pAngVel0(1,:),'ydata',pAngVel0(2,:),'zdata',pAngVel0(3,:));
        set(hAngMom,'xdata',pAngMom0(1,:),'ydata',pAngMom0(2,:),'zdata',pAngMom0(3,:));
        set(hTitle,'string',sprintf('t = %4.2f',t(i)));
        
        if (firsttime)
            firsttime=0;
            pause(0.5);
        end
        drawnow;  
        frame = getframe(gcf); %puts frames in movie
        writeVideo(myV,frame);
    end
end
close(myV); %closes the video formally

function xdot = f(t,x,J1,J2,J3,Jw,n,A,B,x0,Q,R,k)
R = XtoR(x(1:9,1));
w = x(10:end,1);
% %%

% YOUR CODE HERE TO COMPUTE Rdot and wdot 
Rdot = R*skew(w);


w = [w(1); w(2); w(3)]; %state

wdot = (A-B*k)*w;

xdot = [RtoX(Rdot); wdot];

function S = skew(w)
S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];

function X = RtoX(R)
X = reshape(R,9,1);

function R = XtoR(X)
R = reshape(X,3,3);

function [verts,faces,colors] = makebox(x,y,z)
verts = [0 x x 0 0 x x 0; 0 0 0 0 y y y y; 0 0 z z 0 0 z z];
verts = verts - repmat([x/2; y/2; z/2],1,size(verts,2));
faces = [1 2 3 4; 2 6 7 3; 6 5 8 7; 5 1 4 8; 4 3 7 8; 5 6 2 1];
colors(:,:,1) = 1*ones(1,size(faces,1));
colors(:,:,2) = 1*ones(1,size(faces,1));
colors(:,:,3) = 0*ones(1,size(faces,1));
colors(1,5,1) = 1;
colors(1,5,2) = 0;
colors(1,5,3) = 0;
colors(1,3,1) = 0;
colors(1,3,2) = 0;
colors(1,3,3) = 1;
colors(1,2,1) = 0;
colors(1,2,2) = 1;
colors(1,2,3) = 0;




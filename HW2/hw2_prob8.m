function hw2_prob8


% DEFINE THE TIME INTERVAL
tMax = 6;
tRate = 30;
t = linspace(0,tMax,tMax*tRate);

% %%
%
% YOUR CODE HERE TO DEFINE THE INITIAL CONDITIONS
xH0 = zeros(3,1);
R0 = eye(3);
xQ0 = [1 0 0 0];
%
% %%

% INTEGRATE THE ANGULAR RATE EQUATIONS
[t,xR] = ode45(@fR,t,RtoX(R0));
[t,xH] = ode45(@fH,t,xH0);
[t,xQ] = ode45(@fQ,t,xQ0);

figure(3);
subplot(3,1,1);
plot(t,xR);
ylabel('xR')
subplot(3,1,2);
plot(t,xH);
ylabel('xH')
subplot(3,1,3);
plot(t,xQ);
ylabel('xQ')
xlabel('time')
title('Steve Macenski')

% FIND THE RATES
for i=1:length(t)
    rdot(:,i) = RtoX(getRdot(t(i),XtoR(xR(i,:)')));
    thetadot(:,i) = getHdot(t(i),xH(i,:)');
    [q0dot(i),qdot(:,i)] = getQdot(t(i),xQ(i,1),xQ(i,2:4)');
end

% PLOT THE RATES
figure(1);
clf;
subplot(3,1,1);
plot(t,rdot);
ylabel('rdot')
subplot(3,1,2);
plot(t,thetadot);
ylabel('thetadot')
subplot(3,1,3);
plot(t,q0dot,t,qdot);
ylabel('quaterion')
xlabel('time')
title('Steve Macenski')


% CREATE A BOX
[p1,faces,colors] = makebox(0.5,1,0.2);

% SETUP THE PLOT
figure(2);
clf;
axis(1.25*[-1 1 -1 1 -1 1]);
axis equal;
hold on;
plot3(0,0,0,'k.','markersize',16);
hBox = patch('Vertices',p1','Faces',faces,...
          'CData',colors,'FaceColor','flat');
hTitle = title(sprintf('t = %4.2f',0));
lighting flat
light('Position',[0 -2 -1])
light('Position',[0 -2 1])

% ANIMATE THE RESULTS
x = xR;
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
        set(hTitle,'string',sprintf('t = %4.2f',t(i)));
        
        if (firsttime)
            firsttime=0;
            pause(0.5);
        end
        drawnow;
    end
end


function R=HtoR(H)
% %%
%
% YOUR CODE HERE TO CONVERT FROM XYZ EULER ANGLES TO A ROTATION MATRIX
c1 = cos(H(1));
c2 = cos(H(2));
c3 = cos(H(3));
s1 = sin(H(1));
s2 = sin(H(2));
s3 = sin(H(3));


Rxyz = [c2*c3          -c2*s3         s2; ...
        c1*s3+c3*s1*s2 c1*c3-s1*s2*s3 -c2*s1;...
        s1*s3-c1*c3*s2 c3*s1+c1*s2*s3 c1*c2];

R=Rxyz;
%
% %%


function R=QtoR(q0,q)
% %%
%
q = [qo q];
[theta1 theta2 theta3] = quat2angle(q,'XYZ')

c1 = cos(theta1);
c2 = cos(theta2);
c3 = cos(theta3);
s1 = sin(theta1);
s2 = sin(theat2);
s3 = sin(theta3);

Rxyz = [c2*c3          -c2*s3         s2; ...
        c1*s3+c3*s1*s2 c1*c3-s1*s2*s3 -c2*s1;...
        s1*s3-c1*c3*s2 c3*s1+c1*s2*s3 c1*c2];

R=Rxyz;


function [q0dot, qdot] = getQdot(t,q0,q)
% %%
%
% YOUR CODE HERE TO COMPUTE Qdot = (q0dot,qdot)
%
q0
q
w = 10*exp(-t)*[sin(t); sin(2*t); sin(3*t)];


q0dot = -1/2*w'*q
qdot = 1/2.*(w*q0 - [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0]*q)

%
% %%

function thetadot = getHdot(t,theta)
% %%
%
w = 10*exp(-t)*[sin(t); sin(2*t); sin(3*t)];
c1 = cos(theta(1));
c2 = cos(theta(2));
c3 = cos(theta(3));
s1 = sin(theta(1));
s2 = sin(theta(2));
s3 = sin(theta(3));

B = [c3/c2 -s3/c2 0; s3 c3 0; -s2*c3/c2 s2*s3/c2 1];
thetadot = B*w;
%
% %%

function Rdot = getRdot(t,R)
% %%
%
% YOUR CODE HERE TO COMPUTE Rdot
%
w = 10*exp(-t)*[sin(t); sin(2*t); sin(3*t)];

Rdot = R*[0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
%
% %%

function xdot = fQ(t,x)
q0 = x(1,1);
q = x(2:4,1);
[q0dot, qdot] = getQdot(t,q0,q);
xdot = [q0dot; qdot];

function xdot = fH(t,x)
xdot = getHdot(t,x);

function xdot = fR(t,x)
R = XtoR(x);
xdot = RtoX(getRdot(t,R));

function X = RtoX(R)
X = reshape(R,9,1);

function R = XtoR(X)
R = reshape(X,3,3);

function [verts,faces,colors] = makebox(x,y,z)
verts = [0 x x 0 0 x x 0; 0 0 0 0 y y y y; 0 0 z z 0 0 z z];
faces = [1 2 3 4; 2 6 7 3; 6 5 8 7; 5 1 4 8; 4 3 7 8; 5 6 2 1];
colors(:,:,1) = 1*ones(1,size(faces,1));
colors(:,:,2) = 1*ones(1,size(faces,1));
colors(:,:,3) = 0*ones(1,size(faces,1));
colors(1,5,1) = 1;
colors(1,5,2) = 0;
colors(1,5,3) = 0;


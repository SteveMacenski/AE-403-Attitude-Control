
function [R quat]=angle_to_rot
% %%
%
theta = [1.06570000000000;-0.639600000000000;4.25680000000000]
c1 = cos(theta(1));
c2 = cos(theta(2));
c3 = cos(theta(3));
s1 = sin(theta(1));
s2 = sin(theta(2));
s3 = sin(theta(3));

Rxyz = [c2*c3          -c2*s3         s2; ...
        c1*s3+c3*s1*s2 c1*c3-s1*s2*s3 -c2*s1;...
        s1*s3-c1*c3*s2 c3*s1+c1*s2*s3 c1*c2];

R=Rxyz;
quat = angle2quat(1.0657,-.6396,4.2568,'XYZ');

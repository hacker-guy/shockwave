clear all

B= linspace(0,pi,100);
%M = 5, theta = 20
theta1 = 20*(pi/180);

f1 = (2*cot(B)*( ((5^2)*(sin(B).^2)-1)/((5^2)*(1.4+cos(2*B))+2) ) == ...
    tan(theta1));

%M = 5, theta = 30
theta2 = 30*(pi/180);

f2 = 2*cot(B)*( ((5^2)*(sin(B).^2)-1)/((5^2)*(1.4+cos(2*B))+2) ) - tan(theta2);

%M = 5, theta = 45
theta3 = 45*(pi/180);

f3 = 2*cot(B)*( ((5^2)*(sin(B).^2)-1)/((5^2)*(1.4+cos(2*B))+2) ) - ...
    tan(theta3);

plot(B,f2)
ylim([-3 3])
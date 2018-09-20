clear all

syms B

%M = 1.5, theta = 5
theta1 = 5*(pi/180);
B1=0:0.1:pi;

f1 = (2*cot(B1)*( (((1.5)^2)*(power(sin(B1),2))-1)/((1.5^2)*...
    (1.4+cos(2*B1))+2))) - tan(theta1);

%solve(2*cot(B)*( (((M)^2)*(sin(B).^2)-1)/((M^2)*(1.4+cos(2*B))+2))...
%- tan(theta), B)

%M = 1.5, theta = 10
theta2 = 10*(pi/180);
B2=0:0.1:pi;

f2 = (2*cot(B2)*( (((1.5)^2)*(power(sin(B2),2))-1)/((1.5^2)*...
(1.4+cos(2*B2))+2))) - tan(theta2);

%M = 1.5, theta = 15
theta3 = 15*(pi/180);
B3=0:0.1:pi;

f3 = (2*cot(B3)*( (((1.5)^2)*(power(sin(B3),2))-1)/((1.5^2)*...
    (1.4+cos(2*B3))+2))) - tan(theta3);

plot(B1,f1, B2, f2, B3, f3)

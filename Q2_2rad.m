clear all

%M = 1.5, theta = 5
theta1 = 5*(pi/180);
B1=theta1:0.01:(pi/2);

f1 = 2*cot(B1)*( ((1.5^2)*(sin(B1).^2)-1)/((1.5^2)*(1.4+cos(2*B1))+2) ) - ...
    tan(theta1);

%M = 1.5, theta = 10
theta2 = 10*(pi/180);
B2=theta2:0.01:(pi/2);

f2 = 2*cot(B2)*( ((1.5^2)*(sin(B2).^2)-1)/((1.5^2)*(1.4+cos(2*B2))+2) ) - ...
    tan(theta2);

%M = 1.5, theta = 15
theta3 = 15*(pi/180);
B3=theta3:0.01:(pi/2);

f3 = 2*cot(B3)*( ((1.5^2)*(sin(B3).^2)-1)/((1.5^2)*(1.4+cos(2*B3))+2) ) - ...
    tan(theta3);

plot(B1,f1, B2, f2, B3, f3)
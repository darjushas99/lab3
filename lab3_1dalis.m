clear all;
close all
clc;
x = [0.1:1/22:1];
d = (1+0.6*sin(2*pi*x/0.7)+0.3*sin(2*pi*x))/2;

b = randn(1);
w1 = randn(1);
w2 = randn(1);
w = [w1 w2];
niu = 0.05;

c = [0.19 0.92]; % 0.17 0.92
r = [0.18 0.2];  % 0.15 0.19

e = 100;
while abs(e)~=0
   for n=1:length(x)
    y = gauss(x(n),r(1),c(1))*w(1)+gauss(x(n),r(2),c(2))*w(2)+b;
    er = d(n) - y; 
    b = b + niu*er;
    w(1) = w(1) + niu*er*gauss(x(n),r(1),c(1));
    w(2) = w(2) + niu*er*gauss(x(n),r(2),c(2));
   end
   e = e - 1;
end

yy = gauss(x,r(1),c(1))*w(1)+gauss(x,r(2),c(2))*w(2)+b;
plot(x,d,'b',x,yy,'*'); grid on;
legend('Pries apmokyma','Po apmokymo');
 
function F = gauss(x,r,c)
F = exp(-(x-c).^2/(2*r.^2));
end

 
 

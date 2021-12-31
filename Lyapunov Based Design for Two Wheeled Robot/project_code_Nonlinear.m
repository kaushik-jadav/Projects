init = 0.01;
last = 50;
N = 100;
tf =10;
x1_IC = init + (last+init)*randn(N,1);
x2_IC = init + (last+init)*randn(N,1);
theta_d = 10*init + (last+init)*randn(N,1);
%theta_d = 10;
m = 0.1;
g = 9.8;
l = 20;
I = 13.04;
k1 = (m*g*l)/(2*I);
k2 = 1/I;
beta_p = 10;
alpha_p =10;
%b = init + (last+init)*rand(N,1);
%c = init + (last+init)*rand(N,1);
out = sim('proj_ref')
figure(1)
subplot(2,1,1)
plot(out.e)
hold on
xlabel('time');
ylabel('Tracking Error (e)');
grid on
title('Tracking Error')
subplot(2,1,2)
plot(out.r)
hold on
xlabel('time');
ylabel('Filtered Tracking Error (r)');
grid on
title('Filtered Tracking Error') 

figure(2)
subplot(2,1,1)
plot(out.x1)
hold on
xlabel('time');
ylabel('Pitching Angle');
grid on
title('Pitching Angle')

subplot(2,1,2)
plot(out.x2)
hold on
xlabel('time');
ylabel('Pitching Rate');
grid on
title('Pitching Rate')

figure(3)
plot(out.torque)
hold on
xlabel('time');
ylabel('torque');
title('Torque')
grid on


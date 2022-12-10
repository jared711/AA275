clc
clear all
close all

setenceladusglobal
load('/media/jared711/SeagatePortableDrive/UBUNTU/StanfordMATLAB/AA275_Navigation/project/data/NRHO200_eig_idx6_alpha1e-03_QPO.mat')
load("states_trial5.mat")

figure()
plot_rv(xx)
plot_sec

for i = 1:99
    [tt,xx] = ode78e(@(t,x) CR3BP(t,x), 0, T, u(:,i) ,eps);
    plot3(xx(:,1),xx(:,2),xx(:,3),'Color',[.7,.7,.7],'handlevisibility', 'off')
end

Sigma = 1e-9*eye(6);
R = chol(Sigma);
xx_noisy = xx + randn(length(tt),6)*R;

plot_rv(xx_noisy, 'r-x')
legend("True")
legend("Estimated")


%% Plotting altitude
figure()
alt_predict = vecnorm(state_history(1:3,:) - [1-mu;0;0])*RUNIT - SEC.radius;
alt_true = vecnorm(x(1:3,:) - [1-mu;0;0])*RUNIT - SEC.radius;
plot(times*TUNIT/3600,alt_predict,'k')
hold on
plot(times*TUNIT/3600, alt_true,'r')
title('Altitude')
legend('Prediction','Truth','Location','northwest')
xlabel('Time [hrs]')
ylabel('Altitude [km]')

% plotting altitude error
figure()
err = alt_predict - alt_true;
plot(times*TUNIT/3600,err,'k')
xlabel('Time [hrs]')
ylabel('Altitude Error [km]')

title('Altitude Error')


%% Plotting trajectory
figure()
plot_rv(state_history,'r')
plot_rv(x(1:6,:),'k')
plot_sec
view(3)
legend('Prediction','Truth','Location','northeast')

%% Plotting error
figure()
err = vecnorm(state_history(1:3,:) - x(1:3,1:length(state_history)));
plot(times(1:length(state_history))*TUNIT/3600,err*RUNIT,'k')
xlabel('Time [hrs]')
ylabel('Position Error [km]')
title('Position Error')

figure()
err = vecnorm(state_history(4:6,:) - x(4:6,1:length(state_history)));
plot(times(1:length(state_history))*TUNIT/3600,err*VUNIT*1000,'k')
xlabel('Time [hrs]')
ylabel('Velocity Error [m/s]')
title('Velocity Error')

%% Periapsis times
T/2*TUNIT/3600:T*TUNIT/3600:5*T*TUNIT/3600

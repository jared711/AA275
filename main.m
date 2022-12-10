clc
clear all
close all

setenceladusglobal
load('data/NRHO200_eig_idx6_alpha1e-03_QPO.mat')
load("states_trial5.mat")

N_t = 80;

% True Integrated State
% tf = 5*T;
% [tt,xx] = ode78e(@(t,x) CR3BP(t,x), 0, tf, u(:,1) ,eps);


% Attitude (direction cosine matrix)
% q = zeros(q_size,length(tt)); 
% for i = 1:N_t
% %     q(:,i) = get_q(xx(i,1:3)');
%     hs = plot_state(x(:,i));
%     pause(0.1)
%     delete(hs)
% end


% dt = 600/TUNIT; %[NON] 60 sectonds between pictures
% times = 0:dt:tf; % times for pictures to be taken
% N_t = length(times);

% %%%%%%%%% PLACEHOLDER %%%%%%%%%%%
% x_interp = interp1(tt, xx(:,1), times);
% y_interp = interp1(tt, xx(:,2), times);
% z_interp = interp1(tt, xx(:,3), times);
% r_interp = [x_interp;
%             y_interp;
%             z_interp];
% q_interp = zeros(q_size, N_t);
% for i = 1:q_size
%     q_interp(i,:) = interp1(tt, q(i,:), times);
% end
q_size = 0;

view_mat_corr = [1,0,0;0,0,1;0,-1,0];

global r_err v_err
r_err = 0.1/RUNIT;
v_err = 0.0001/VUNIT;
% initial prediction taken from a Gaussian distribtution around truth
Q = blkdiag(0.1/RUNIT*eye(3), 0.01/VUNIT*eye(3), eps*eye(q_size));
R = blkdiag(0.01/RUNIT*eye(3), eps*eye(q_size)); % process noise
n = length(Q); % state size 6 for r+v, 9 for direction cosine matrix
m = length(R); % measurement size

% initial estimated state
% x_t = [xx(1,:)'; q(:,1)] + chol(Q)'*randn(n,1);
x_t = x(1:6,1);% + chol(Q)'*randn(n,1);

state_history = zeros(6,N_t);
state_history(:,1) = x_t(1:6);

% Kalman Filter
P_t = Q;
figure()
plot_rv(xx,'k',4,0.5,'off')
hold on
plot_sec

% Prepare the new file.
% vidObj = VideoWriter('altonly');
% open(vidObj);

for t = 1:N_t-1
    fprintf("t = %i\n",t)

    hs = plot_state(x_t, x(:,t));
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);

    pause(0.1)
    delete(hs)
    x_t1 = x_t; % Current prediction becomes old prediction
    P_t1 = P_t;

    x_t_t1 = f(x_t1, dt);
    F_t = dfdx(x_t1, dt); % Jacobian of state transition function
    P_t_t1 = F_t*P_t1*F_t + Q;

    if norm(x(1:3,t+1)-[1-mu;0;0])*RUNIT < 500
        z_t = get_state(x(:,t+1));
    else
        z_t = get_vector_from_camera(t+1, view_mat(:,:,t+1)); % actual measurement
    end
    % Let's make a different function for measurement and preduicted
    % measurement
    H_t = dhdx(x_t1, m, n); % Jacobian of measurement function
    K_t = P_t_t1*H_t'*inv(R + H_t*P_t_t1*H_t');
    x_t = x_t_t1 + K_t*(z_t - h(x_t_t1));
    P_t = (eye(n) - K_t*H_t)*P_t_t1;

    state_history(:,t+1) = x_t(1:6);
end

% Close the file.
% close(vidObj);

% predicted measurement
function z = h(x)
global mu RUNIT
z = x(1:3);
end

% predicted state
function x_t_t1 = f(x,dt)
[~,xx] = ode78e(@(t,x) CR3BP(t,x), 0, dt, x(1:6) ,eps);
rv = xx(end,:)';
% q = get_q(rv(1:3));
% x_t_t1 = [rv;q];
x_t_t1 = rv;
end

function H_t = dhdx(x, m, n)
global mu RUNIT
% if norm(x(1:3)-[1-mu;0;0])*RUNIT < 500
%     re = [1-mu;0;0];   
%     r = x(1:3);
%     H_t = zeros(1,n);
%     H_t(1,1:3) = unit(r-re)';
% else
H_t = [eye(3), zeros(3)];
% end
% H_t(2:m, 7:n) = eye(n-6);
end

function F_t = dfdx(x, dt)
n = length(x);
rv = x(1:6);
F_t = zeros(n);
% F_t(1:6,1:6) = CR3BP_Jacobian(rv);
s = [reshape(eye(6),36,1); rv];
[tt,xx] = ode78e(@(t,x) CR3BP_STM(t,x), 0, dt, s ,eps);
F_t2(1:6,1:6) = reshape(xx(end,1:36),6,6);
end

function q = get_q(r)
global mu
uhat = unit([1-mu;0;0] - r);
vhat = unit(cross(uhat, [1;0;0]));
what = unit(cross(uhat, vhat));
A = [uhat, vhat, what]; % A transforms from spacecraft frame to Enceladus rotating frame
q = reshape(A, numel(A), 1);
end

function hs = plot_state(x_predict, x_true)
h1 = plot_rv(x_predict(1:3),'rx');
h2 = plot_rv(x_true(1:3),'kx');
legend("Prediction","Truth")
view(3)
hs = [h1,h2];
% h2 = plot_quiver(x_t(1:6));
% h3 = plot_quiver([x_t(1:3);0.01*x_t(7:9)]);
% h4 = plot_quiver([x_t(1:3);0.01*x_t(10:12)]);
% h5 = plot_quiver([x_t(1:3);0.01*x_t(13:15)]);
% hs = [h1,h2,h3,h4,h5];
end

function alt = get_alt(x_t)
global mu SEC RUNIT r_err
re = [1-mu;0;0];
alt = norm(x_t(1:3) - re) - SEC.radius/RUNIT;
alt = alt + r_err/RUNIT*randn;
end

function p = get_point(x_t)
global mu RUNIT r_err
re = [1-mu;0;0];
p = unit(re - x_t(1:3));
p = p + r_err/RUNIT*randn(3,1);
end

function s = get_state(x_t)
global SEC mu RUNIT
alt = get_alt(x_t);
p = get_point(x_t);
s = -(alt+SEC.radius/RUNIT)*p + [1-mu;0;0];
end

function alt = altimeter(x_t)
global mu SEC RUNIT
re = [1-mu;0;0];
alt = norm(x_t(1:3) - re) - SEC.radius/RUNIT; % [NON]
sigma = 0.05;% 5% error on altitude measurements from altimeter
alt = alt*(1 + sigma*randn); 
end
    
function [C_mat, p] = get_ellipse(t)
global SEC RUNIT
% x_t is the true state with attitude in it
I = imread(sprintf("fig/trial5/img_%.04i.jpg",t));
Igray = rgb2gray(I);
[y,x] = find(Igray);
a = (max(x)-min(x))/2;
b = (max(y)-min(y))/2;
% xc = (max(x)+min(x))/2;
% yc = (max(y)+min(y))/2;
% r = x_t(1:3)
A = 1/a^2;
B = 0;
C = 1/b^2;
D = 0;
F = 0;
G = A;
f = 7.153492551587466e+02;
C_mat = [  A*f^2, B*f^2/2, D*f/2;
         B*f^2/2,   C*f^2, F*f/2;
           D*f/2,   F*f/2,     G];
p = f*SEC.radius/RUNIT/a;

end

function r = get_vector_from_camera(t, view_mat)
global SEC RUNIT mu
T_C2P = view_mat(1:3,1:3); % rotation matrix from camera frame to planet frame
ellipse_axes = SEC.radius/RUNIT*eye(3);
Ainv = T_C2P' * ellipse_axes * T_C2P;
[C,p] = get_ellipse(t);
A = inv(Ainv);

[V,Lambda] = eig(Ainv*C);
signs = sign(Lambda);
if signs(1) == signs(2)
    idx = 3;
elseif signs(1) == signs(3)
    idx = 2;
elseif signs(2) == signs(3)
    idx = 1;
end
lambda = Lambda(idx);
v = V(:,idx);

rho2 = (trace(C) - lambda*trace(A))/(lambda*v'*A*A*v - lambda*v'*A*v*trace(A));
r = -sign(v(3)) * sqrt(rho2)*T_C2P*v;
r = -T_C2P'*[0;0;p] + [1-mu;0;0];
end


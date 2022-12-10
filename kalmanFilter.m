
global mu

xbar_t = f(x_t1) + w_t;
z_t = h(x_t) + v_t;

% % % % x_t = % true state
% % % % xbar_t % estimated state
% % % % z_t % measurement

for t = 1:N_t


    x_t_t1 = f(x_t1);
    F_t = dfdx(x_t1);
    P_t_t1 = F_t*P_t1*F_t + Q;

    H_t = dhdx(x_t1);

    K_t = P_t_t1*H_t'*inv(R + H_t*P_t_t1*H_t);
    x_t = x_t_t1 + K_t*(z_t - h(x_t_t1));
    P_t = (eye(10) - K_t*H_t)*P_t_t1;

end

function z = h(x)
re = [1-mu;0;0];   
r = x(1:3);
d = norm(r-re);
q = x(7:10);
z = [d;q];
end

function x_t_t1 = f(x,dt)
[~,xx] = ode78e(@(t,x) CR3BP(t,x), 0, dt, x ,eps);
x_t_t1 = xx(end,:)';
end

function H_t = dhdx(x)
re = [1-mu;0;0];   
r = x(1:3);
H_t = [unit(r-re)', zeros(7,1)];
end

function F_t = dfdx(x)
rv = x(1:6);
F_t = zeros(10);
F_t(1:6,1:6) = CR3BPJac(rv);
end
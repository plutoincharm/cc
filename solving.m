clearvars; clc;
global g k1 k2 l1 l2 l3 l4 m1 m2 m3 m4 M D c1 c2;    

% constants required:
g = 9.81;
k1 = 5;
k2 = 5;
l1 = 0.4418;
l2 = 0.4033;

r1 = 0.1872;
r2 = 0.1738;


% l3 = 1;
% l4 = 1;
D = 60;
c1 = 2;
c2 = 2;

m = [7.4 3.41 1.073/2 1.073/2];
m1 = m(1);
m2 = m(2);
m3 = m(3);
m4 = m(4);
M = 50.172;

% parameters
tf = 10;
tspan = [0,tf];
dt = 1e-3;
t = 0:dt:tf;
n = length(t);
% variable order: (omg1, omg2, omg3, xdot, ydot, tht1, tht2, tht3, x, y)
q0 = [0;0;0;0;0;0;0;30;0; (m1 + m2 + m3 + m4) * g / k2];

handle_f = @pendulum_dynamics;  

for i = 1:n-1
 if i==1
 q(i,:) = q0';
 qnew = fun_explicitEuler(handle_f,q(i,:)',dt,t(i));
 q(i+1,:) = qnew';
 else
 qnew = fun_explicitEuler(handle_f,q(i,:)',dt,t(i));
 q(i+1,:) = qnew';
 end
end

qExEu = q;
% save('SlidePenduLagODE_ExEu.mat','t','q');

tiledlayout(3,4);

for i = 1:10
    nexttile
    plot(t,qExEu(:,i));
    names = (["omg1", "omg2", "omg3", "xdot", "ydot", "tht1", "tht2", "tht3", "x", "y"]);
    title(strcat('Explicit Euler:', names(i)));
    xlabel('t \rightarrow');
end

function val = T1(t)
    val = 0;
end

function val = T2(t)
    val = 0;
end

function val = T3(t)
if (t < 2) || (t > 3)    
    val = 0;
else
    val = 1;
end
end


function qdot = pendulum_dynamics(q,t)
    global g k1 k2 l1 l2 l3 l4 m1 m2 m3 m4 M D c1 c2;    
    
  
    omg1 = q(1);
    omg2 = q(2);
    omg3 = q(3);
    xdot = q(4);
    ydot = q(5);
    tht1 = q(6);
    tht2 = q(7);
    tht3 = q(8);
    x = q(9);
    y = q(10);
    
    A1 = [1/12 * l1^2 * (7 * m1 + 12 * (m2 + m3 + m4)), 1/2 * l1 * l2 * (m2 + 2 * (m3 + m4)) * cos(tht1 - tht2), 1/2 * l1 * (l3 * m3 * cos(tht1 - tht3) + l4 * m4 * cos(D - tht1 + tht3)), 1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * cos(tht1), 1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * sin(tht1)];
    A2 = [1/2 * l1 * l2 * (m2 + 2 * (m3 + m4)) * cos(tht1 - tht2), 1/12 * l2^2 * (7 * m2 + 12 * (m3 + m4)), 1/2 * l2 * (l3 * m3 * cos(tht2 - tht3) + l4 * m4 * cos(D - tht2 + tht3)), 1/2 * l2 * (m2 + 2 * (m3 + m4)) * cos(tht2), 1/2 * l2 * (m2 + 2 * (m3 + m4)) * sin(tht2)];
    A3 = [1/2 * l1 * (l3 * m3 * cos(tht1 - tht3) + l4 * m4 * cos(D - tht1 + tht3)), 1/2 * l2 * (l3 * m3 * cos(tht2 - tht3) + l4 * m4 * cos(D - tht2 + tht3)), (7 * l3^2 * m3)/12 + (7 * l4^2 * m4)/12, 1/2 * l3 * m3 * cos(tht3) + 1/2 * l4 * m4 * cos(D + tht3), 1/2 * l3 * m3 * sin(tht3) + 1/2 * l4 * m4 * sin(D + tht3)];
    A4 = [1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * cos(tht1), 1/2 * l2 * (m2 + 2 * (m3 + m4)) * cos(tht2), 1/2 * l3 * m3 * cos(tht3) + 1/2 * l4 * m4 * cos(D + tht3), M + m1 + m2 + m3 + m4, 0];
    A5 = [1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * sin(tht1), 1/2 * l2 * (m2 + 2 * (m3 + m4)) * sin(tht2), 1/2 * l3 * m3 * sin(tht3) + 1/2 * l4 * m4 * sin(D + tht3), 0, M + m1 + m2 + m3 + m4];
        
    A = -[A1; A2; A3; A4; A5];
    
    Mmat = [A,zeros(5);zeros(5),eye(5)];
    
    b1 = -T1(t) + T2(t) + T3(t) + 1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * omg1 * ydot * cos(tht1) + 1/2 * g * l1 * (m1 + 2 * (m2 + m3 + m4)) * sin(tht1) - 1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * omg1 * xdot * sin(tht1) - 1/2 * l1 * l2 * (m2 + 2 * (m3 + m4)) * (omg1 - omg2) * omg2 * sin(tht1 - tht2) + 1/2 * l1 * omg3 * (-omg1 + omg3) * (l3 * m3 * sin(tht1 - tht3) - l4 * m4 * sin(D - tht1 + tht3)) - 1/2 * l1 * omg1 * ((m1 + 2 * (m2 + m3 + m4)) * (ydot * cos(tht1) - xdot * sin(tht1)) - l2 * (m2 + 2 * (m3 + m4)) * omg2 * sin(tht1 - tht2) + omg3 * (-l3 * m3 * sin(tht1 - tht3) + l4 * m4 * sin(D - tht1 + tht3)));
    b2 = -T2(t) + T3(t) + 1/2 * l2 * (m2 + 2 * (m3 + m4)) * omg2 * ydot * cos(tht2) - 1/2 * l1 * l2 * (m2 + 2 * (m3 + m4)) * omg1 * (omg1 - omg2) * sin (tht1 - tht2) + 1/2 * g * l2 * (m2 + 2 * (m3 + m4)) * sin(tht2) - 1/2 * l2 * (m2 + 2 * (m3 + m4)) * omg2 * xdot * sin(tht2) + 1/2 * l2 * omg3 * (-omg2 + omg3) * (l3 * m3 * sin(tht2 - tht3) - l4 * m4 * sin(D - tht2 + tht3)) - 1/2 * l2 * omg2 * (l1 * (m2 + 2 * (m3 + m4)) * omg1 * sin(tht1 - tht2) + (m2 + 2 * (m3 + m4)) * (ydot * cos(tht2) - xdot * sin(tht2)) + omg3 * (-l3 * m3 * sin(tht2 - tht3) + l4 * m4 * sin(D - tht2 + tht3))); 
    b3 = 1/2 * l3 * m3 * omg3 * ydot * cos(tht3) + 1/2 * l4 * m4 * omg3 * ydot * cos(D + tht3) - 1/2 * l3 * m3 * omg3 * xdot * sin(tht3) - 1/2 * l4 * m4 * omg3 * xdot * sin(D + tht3) + 1/2 * l1 * omg1 * (omg1 - omg3) * (-l3 * m3 * sin(tht1 - tht3) + l4 * m4 * sin(D - tht1 + tht3)) + 1/2 * l2 * omg2 * (omg2 - omg3) * (-l3 * m3 * sin(tht2 - tht3) + l4 * m4 * sin(D - tht2 + tht3)) + 1/2 * (-2 * T3(t) + g * l3 * m3 * sin(tht3) + g * l4 * m4 * sin(D + tht3) - omg3 * (l3 * m3 * ydot * cos(tht3) + l4 * m4 * ydot * cos(D + tht3) - l3 * m3 * xdot * sin(tht3) - l4 * m4 * xdot * sin(D + tht3) + l1 * omg1 * (l3 * m3 * sin(tht1 - tht3) - l4 * m4 * sin(D - tht1 + tht3)) + l2 * omg2 * (l3 * m3 * sin(tht2 - tht3) - l4 * m4 * sin(D - tht2 + tht3)))); 
    b4 = k1 * x + c1 * xdot + 1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * omg1^2 * sin(tht1) - 1/2 * l2 * (m2 + 2 * (m3 + m4)) * omg2^2 * sin(tht2) - 1/2 * l3 * m3 * omg3^2 * sin(tht3) - 1/2 * l4 * m4 * omg3^2 * sin(D + tht3);
    b5 = -g * (m1 + m2 + m3 + m4) + k2 * y + c2 * ydot + 1/2 * l1 * (m1 + 2 * (m2 + m3 + m4)) * omg1^2 * cos(tht1) + 1/2 * l2 * (m2 + 2 * (m3 + m4)) * omg2^2 * cos(tht2) + 1/2 * l3 * m3 * omg3^2 * cos(tht3) + 1/2 * l4 * m4 * omg3^2 * cos(D + tht3);
    
    b = [b1; b2; b3; b4; b5; omg1; omg2; omg3; xdot; ydot];
    
    qdot = Mmat\b;

end

function y_n_plus1 = fun_explicitEuler(handle_f,y_n,dt,t)
    global g k1 k2 l1 l2 l3 l4 m1 m2 m3 m4 M D c1 c2;    
 
   
    if nargin==3 
     y_n_plus1 = y_n + dt*feval(handle_f,y_n);
    else
     y_n_plus1 = y_n + dt*feval(handle_f,y_n,t);
    end
end


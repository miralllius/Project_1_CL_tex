% Project 1_2_a, plot solution for (5) = 2, (6) = 3
close all
clear all
clc

%% Initialize

IC = 3;

num_flux_1 = 'LF'; % for project 1_2_a
num_flux_2 = 'Roe'; % for project 1_3_c

[U0, S, a, b, bc, g] = Initial_conditions(IC);
CFL = 0.5; T = 2;
N = 1000;

%Load reference solution
switch IC
    case 1
        h = (b-a)/N;
        xf = a:h:b;
        xc_ref = a+0.5*h:h:b-0.5*h;
        U_ref = zeros(2,N);
        for i = 1:N
            U_ref(:,i) = integral(U_ex, xf(i), xf(i+1), 'AbsTol', 1e-14, 'ArrayValued', true)/h;
        end
    case 2
        load('Ref_IC2');
    case 3
        load('Ref_IC3');
    case 4
        load('Ref_IC4');
end

xc_ref = xc;
U_ref = U1;
%% Compute the solution
U1 = solver(U0,S,a,b,N,T,CFL,bc,num_flux_1);
U2 = solver(U0,S,a,b,N,T,CFL,bc,num_flux_2);

%% Plot the solution

h = (b-a)/N;
xc = a+0.5*h:h:b-0.5*h;

figure()
subplot(2,1,1)
plot(xc_ref,U_ref(1,:), '-k', 'linewidth', 2)
hold on
plot(xc, U1(1,:), '-', 'LineWidth', 2)
plot(xc, U2(1,:), '-', 'LineWidth', 2)
ylabel('Height')
xlabel('x')
legend('Reference','LF','Roe','Location','Best')
set(legend,'FontSize',14)

subplot(2,1,2)
plot(xc_ref,U_ref(2,:), '-k', 'linewidth', 2)
hold on
plot(xc, U1(2,:), '-', 'LineWidth', 2)
plot(xc, U2(2,:), '-', 'LineWidth', 2)
ylabel('Discharge')
xlabel('x')
legend('Reference','LF','Roe','Location','Best')
set(legend,'FontSize',14)
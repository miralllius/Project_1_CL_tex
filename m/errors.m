% Project 1_2_b/1_3_c

%% Initialization
IC=1;
num_flux_1 = 'LF'; % for project 1_2_b
num_flux_2 = 'Roe'; % for project 1_3_c

[U0, S, a, b, bc, g] = Initial_conditions(IC);
CFL = 0.5; T = 2;
N_set = [10,20,50,100,500, 1000];

%Reference solution

switch IC
    case 1
        N = 10000;
        U_ex = @(x) U0(x-a*T);
        U_ref = zeros(2,N);
        h = (b-a)/N;
        xf = a:h:b;
        xc_ref = a+0.5*h:h:b-0.5*h;
        for i = 1:N
            display(i)
            U_ref(:,i) = integral(U_ex, xf(i), xf(i+1), 'AbsTol', 1e-14, 'ArrayValued', true)/h;
        end
    case 2
        load('Ref_IC2')
    case 3
        load('Ref_IC3')
        U_ref = U1;
        xc_ref = xc;
    case 4
        load('Ref_IC4')
        T=0.5;
end

errors1 = zeros(2,length(N_set));
errors2 = zeros(2,length(N_set));
p = 2; % to compute p norm

%% Compute the solutions and errors

for i = 1:length(N_set)
    
    N = N_set(i);
    
    U1 = solver(U0,S,a,b,N,T,CFL,bc,num_flux_1);
    U2 = solver(U0,S,a,b,N,T,CFL,bc,num_flux_2);
    
    h = (b-a)/N;
    xc = a+0.5*h:h:b-0.5*h;
    errors1(:,i) = p_error(U1, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors2(:,i) = p_error(U2, ref_to_current(U_ref,xc_ref,xc), h, p);
    
    
end

%% Plot Errors

figure()
%suptitle('Errors')

subplot(2,1,1)
loglog(N_set,N_set.^(-1), '-k', 'LineWidth', 2)
hold on
r1 = polyfit(log(2./N_set),log(errors1(1,:)),1);
loglog(N_set, errors1(1,:), '--', 'LineWidth', 2)
g1 = polyfit(log(2./N_set),log(errors2(1,:)),1);
loglog(N_set, errors2(1,:), '--', 'LineWidth', 2)
grid on
legend(texlabel('log(N^{-1})'), ['LF error, ' sprintf('slope = %f',r1(1)) ],['Roe error, ' sprintf('slope = %f',g1(1))])
set(legend,'FontSize',12)
ylabel('Height error')
xlabel('Number of cells')

subplot(2,1,2)
grid on
loglog(N_set,N_set.^(-1), '-k', 'LineWidth', 2)
hold on
r2 = polyfit(log(2./N_set),log(errors1(2,:)),1);
loglog(N_set, errors1(2,:), '--', 'LineWidth', 2)
g2 = polyfit(log(2./N_set),log(errors2(2,:)),1);
loglog(N_set, errors2(2,:), '--', 'LineWidth', 2)
grid on
legend(texlabel('log(N^{-1})'), ['LF error, ' sprintf('slope = %f',r2(1)) ],['Roe error, ' sprintf('slope = %f',g2(1))])
set(legend,'FontSize',12)
ylabel('Discharge error')
xlabel('Number of cells')




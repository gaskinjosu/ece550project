
%Global variables
delta = (1/20);
RC = 1;

set(groot,'defaultLineMarkerSize',10);
set(groot,'defaultLineLineWidth',1.5);

%2(a) - Verify the approximation by comparing the state vector x(t) and the
%output y(t) (for t0 = 0 <= t <= 20[sec]) to their theoretical values in the 
%case of n = 1 and n = 2 (i.e., 1 ×1 and 2 ×2 state matrix).
%Consider the zero-state case with u(t) = mu(t) (a step function) 
%and the zero-input case with x(0) specified above. We define the step 
%function mu(t) as mu(t) = 1 for t >= 0 and mu(t) = 0 for t < 0.

%IMPORTANT: if we are considering t between 0 and 20, this means that we
%are considering k between 0 and 400 (i.e we are samping 20 times per
%second) -- (use delta*k = t to convert)

%%%%%%%%%%%%% N=1 %%%%%%%%%%%%
A = (1/RC)*[-1];
B = (1/RC)*[1];
C = [1];
D = 0;

A_prime = [1] + A*delta;
B_prime = B*delta;
C_prime = C;
D_prime = D;

%ZERO STATE: If we consider the zero-state resposne,
% then we can see that x[0]= 0 and all other x[k] = A_prime*x[k] + B_prime = 1/RC*delta.
%Similarly, x'(t) = A*x(t) + B*mu(t) = -x(t) + 1 for t>=0 in this setting, which implies that x(t) is
%equal to...integral(0 to t) of e^-t = 1 - e^-t
%Additionally, y[k] = C_prime*x[k] = x[k] and y(t) = C*x(t) = x(t)

%Setting up the time axis for the discrete case 
%(so the points will be plotted as going from 0 to 20 instead of based on array indices)
x = 0:delta:20;
[x_zero_state_discrete_1,y_zero_state_discrete_1] = zeroState(A_prime,B_prime,C_prime,D_prime,1);

%Creating the continuous time equivalents
t = 0:0.01:20;
x_zero_state_continuous_1 = 1 - exp(-t);
y_zero_state_continuous_1 = 1 - exp(-t);

%Plotting the discrete and continuous x on the same axes
figure()
plot(x,x_zero_state_discrete_1,'linestyle','none','marker','.');
hold on;
plot(t,x_zero_state_continuous_1);

%Plotting the discrete and continuous y on the same axes
figure()
plot(x,y_zero_state_discrete_1,'linestyle','none','marker','.');
hold on;
plot(t,y_zero_state_continuous_1);


%ZERO INPUT: If we consider zero-input response (with x(0) = 1) then we
%will obtain x[k+1] = A_prime * x[k] => x[1] = A_prime, X[2] = A_prime^2,
%etc...which is repeated multiplication of A_prime = 1+ -1(1/20) = 19/20.
%This results in a number that continues to get smaller with time
%we also obtain, x'(t) = Ax(t) = -x(t), which implies that x(t) =
%(1)*e^-(t), which also continues to get smaller with time.
%Additionally, y[k] = C_prime*x[k] = x[k] (this makes sense since there is
%only a single capacitor so it is both the state and the observation)
%and y(t) = C*x(t) = x(t) 

x = 0:delta:20;
[x_zero_input_discrete_1,y_zero_input_discrete_1] = zeroInput(A_prime,B_prime,C_prime,D_prime,1);

t = 0:0.001:20;
x_zero_input_continuous_1 = exp(-t);
y_zero_input_continuous_1 = exp(-t);

figure()
plot(x,x_zero_input_discrete_1,'linestyle','none','marker','.');
hold on;
plot(t,x_zero_input_continuous_1);

figure()
plot(x,y_zero_input_discrete_1,'linestyle','none','marker','.');
hold on;
plot(t,y_zero_input_continuous_1);


%%%%%%%%%%%%%% N=2 %%%%%%%%%%%%%%%
A = (1/RC)*[-1 1;1 -2];
B = (1/RC)*[0 1];
C = [1 0];
D = 0;

A_prime = [1 0;0 1] + A*delta;
B_prime = B*delta;
C_prime = C;
D_prime = D;

%ZERO STATE:
x = 0:delta:20;
[x_zero_state_discrete_2,y_zero_state_discrete_2] = zeroState(A_prime,B_prime,C_prime,D_prime,2);

t = 0:0.01:20;
%CONTINUOUS TIME EQUIVALENTS:
%This can be solved for by consdering x(t) = integral(0 to t) of (Phi*B*U)
%(since the initial state is zero). Since input is a step function, U can
%be set = 1. So, only the transition matrix need to be determined.

%Doing an eigendecomposition of A allows us to express it using:
U_2 = [(1+sqrt(5))/2 (1-sqrt(5))/2; 1 1];
Lambda_2 = [(sqrt(5)-3)/2 0; 0 (-sqrt(5)-3)/2];
%A = U_2*Lambda_2*inv(U_2)

%Defining symbolic variables to be used in matrix multiplication
syms t_sym
syms sigma_sym

%Now, the transition matrix can be expressed as Phi = U*e^(lambda*t)*U_inv
%Defining the matrix that is e^(lambda*(t-sigma))
e_lambda_t_sigma = [exp((sqrt(5)-3)*(t_sym - sigma_sym)/2) 0; 0 exp((-sqrt(5)-3)*(t_sym - sigma_sym)/2)];
%Defining Phi(t,sigma)
phi_t_sigma = U_2*e_lambda_t_sigma*inv(U_2);

%Now, finding the integrand for the zero state response:
integrand = phi_t_sigma*(B.');

%Integrating to get the continuous time zero state resposne:
x_zero_state_continuous_2 = int(integrand,sigma_sym,0,t_sym);
y_zero_state_continuous_2 = (C)*x_zero_state_continuous_2;

figure()
plot(x,x_zero_state_discrete_2,'linestyle','none','marker','.');
hold on;
%Using fplot instead of plot here due to symbolic variables being used
fplot(x_zero_state_continuous_2,[0 20]);

figure()
plot(x,y_zero_state_discrete_2,'linestyle','none','marker','.');
hold on;
fplot(y_zero_state_continuous_2,[0 20]);


%ZERO INPUT:
x = 0:delta:20;
[x_zero_input_discrete_2,y_zero_input_discrete_2] = zeroInput(A_prime,B_prime,C_prime,D_prime,2);

t = 0:0.01:20;
%CONTINUOUS TIME EQUIVALENTS:

%We can simply re-use the transition matrix found for the zero-state
%response here, expect now it will be a function of t and t_0 = 0.
e_lambda_t = [exp((sqrt(5)-3)*(t_sym)/2) 0; 0 exp((-sqrt(5)-3)*(t_sym)/2)];
phi_t = U_2*e_lambda_t*inv(U_2);

%Now, the zero-input response is simply the initial state multiplied by
%this transition matrix
x_zero_input_continuous_2 = phi_t*([1 0].') ;
y_zero_input_continuous_2 = (C)*x_zero_input_continuous_2;

figure()
plot(x,x_zero_input_discrete_2,'linestyle','none','marker','.');
hold on;
fplot(x_zero_input_continuous_2,[0 20]);

figure()
plot(x,y_zero_input_discrete_2,'linestyle','none','marker','.');
hold on;
fplot(y_zero_input_continuous_2,[0 20]);



%2(b)
%N=20:
A = (1/RC)*[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0  1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2;];
       
B = (1/RC)*[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
C = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
D = 0;

A_prime = eye(20) + A*delta;
B_prime = B*delta;
C_prime = C;
D_prime = D;

%ZERO STATE:
[x_zero_state_discrete_20,y_zero_state_discrete_20] = zeroState(A_prime,B_prime,C_prime,D_prime,20);

figure()
plot(x,x_zero_state_discrete_20);

figure()
plot(x,y_zero_state_discrete_20);

%ZERO INPUT:
[x_zero_input_discrete_20,y_zero_input_discrete_20] = zeroInput(A_prime,B_prime,C_prime,D_prime,20);

figure()
plot(x,x_zero_input_discrete_20);

figure()
plot(x,y_zero_input_discrete_20);


%4(b): The Lyapunov stability criteria with Q = identity is as follows:
% US : A^T + A <= 0, which means -A^T - A >= 0 (-A^T - A is PSD)
% UES: A^T + A <= -vI which means -vI - A^T - A >= 0 (-vI - A^T - A is PSD)

%Calculating US lyapunov condition
US_lyapunov = (-1*(A.')) + (-1*A);
%printing out eigenvalues, we can see they are all >= 0
eig(US_lyapunov)

%Calculating UES lyapunov condition
%(setting v=0.01 achieves all eigenvalues >= 0)
UES_lyapunov = (-0.01*eye(20)) + (-1*(A.')) + (-1*A);
%printing out eigenvalues, we can see they are all >= 0
eig(UES_lyapunov)


%4(c): for BIBO stability, we can use the argument that if the norm of C
%and the norm of B are bounded (which they obviously are in our case), then
%if the system is UES it is also BIBO stable.
%=> since the system is UES (by Lyapunov UES condition), then the system is
%also BIBO stable


%TODO: for #5 and #6, need to determine if we should be using the discrete
%A_prime, B_prime, C_prime, or the continuous A,B,C for calculation (it does not seem to make
%a difference). Additionally, need to determine if setting the tolerance in
%the rank() function so low is not providing a misleading result (because
%if the tolerance is raised, rank() reports less than rank n)

%5: For n=20, is the system controllable on [0, 20]
%Since the system is LTI, we can compute the controllability matrix
%and ensure that it is rank n to determine if the system is controllable
%Constructing the controllability matrix
current_col = (B.');
controllability_matrix = [current_col];
for i = 2:1:20
    current_col = A*current_col;
    controllability_matrix = [controllability_matrix current_col];
end
%Printing out the rank of the controllability matrix
%note: had to set tolerance to be 1e-8 or smaller to get
%result that controllability matrix is full rank
controllability_rank = rank(controllability_matrix,1e-50);
disp('Controllability Matrix Rank:')
disp(controllability_rank)

%6: For n=20, is the system observable on [0, 20]
%Since the system is LTI, we can compute the observability matrix
%and ensure that it is rank n to determine if the system is observable
%Constructing the observability matrix
current_row = C;
observability_matrix = [current_row];
for i=2:1:20
    current_row = current_row*A;
    observability_matrix = [observability_matrix ; current_row];
end
%Printing out the rank of the controllability matrix
%note: had to set tolerance to be 1e-8 or smaller to get
%result that controllability matrix is full rank
observability_rank = rank(observability_matrix,1e-30);
disp('Observability Matrix Rank:')
disp(observability_rank)



function [x,y] = zeroState(A_prime,B_prime,C_prime,D_prime,n)
    x(:,1) = zeros(1,n);
    %Note: matlab is not zero-indexed, so the first element must be index 1
    for k = 1:1:400
        x(:,k+1) = A_prime*x(:,k) + (B_prime*1).';
        y(:,k) = C_prime*x(:,k) + D_prime*1;
    end
    y(401) = C_prime*x(:,401) + D_prime*1;
end

function [x,y] = zeroInput(A_prime,B_prime,C_prime,D_prime,n)
    x(:,1) = zeros(1,n);
    x(1,1) = 1;
    %Note: matlab is not zero-indexed, so the first element must be index 1
    for k = 1:1:400
        x(:,k+1) = A_prime*x(:,k);
        y(:,k) = C_prime*x(:,k);
    end
    y(401) = C_prime*x(:,401);
end


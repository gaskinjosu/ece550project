
%Global variables
delta = (1/20);
RC = 1;

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
plot(x,x_zero_state_discrete_1,'linestyle','none','marker','o');
hold on;
plot(t,x_zero_state_continuous_1);

%Plotting the discrete and continuous y on the same axes
figure()
plot(x,y_zero_state_discrete_1,'linestyle','none','marker','o');
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
plot(x,x_zero_input_discrete_1,'linestyle','none','marker','o');
hold on;
plot(t,x_zero_input_continuous_1);

figure()
plot(x,y_zero_input_discrete_1,'linestyle','none','marker','o');
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
%TODO: SOLVE ANALYITCALLY AND FILL THIS IN???
%x_zero_state_continuous_2 = ;
%y_zero_state_continuous_2 = ;

figure()
plot(x,x_zero_state_discrete_2,'linestyle','none','marker','o');
%hold on;
%plot(t,x_zero_state_continuous_2);

figure()
plot(x,y_zero_state_discrete_2,'linestyle','none','marker','o');
%hold on;
%plot(t,y_zero_state_continuous_2);


%ZERO INPUT:
x = 0:delta:20;
[x_zero_input_discrete_2,y_zero_input_discrete_2] = zeroInput(A_prime,B_prime,C_prime,D_prime,2);

t = 0:0.01:20;
%TODO: SOLVE ANALYITCALLY AND FILL THIS IN???
%x_zero_input_continuous_2 = ;
%y_zero_input_continuous_2 = ;

figure()
plot(x,x_zero_input_discrete_2,'linestyle','none','marker','o');
%hold on;
%plot(t,x_zero_input_continuous_2);

figure()
plot(x,y_zero_input_discrete_2,'linestyle','none','marker','o');
%hold on;
%plot(t,y_zero_input_continuous_2);



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


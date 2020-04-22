%Script that uses MultifractalV2 with:
%Inputs -     Matriz: Square and 2 exponent size Matrix
%             Q:      Exponents chosen in this script, tipically is Qmax=10

%Outputs -    A:      Matrix of sum(mu*P) with dependance with eps and Q
%             eps:    Array of the sizes of the box, also exponents of 2 
%             tau:    Matrix of the different values of tau with a
%                     dependance in eps and Q
%             IQ:     Matrix of sum(P^q) with dependance in Q and eps (I
%                     dont use it for the rest

% With the ouputs we calculate the slopes and save and show the results. It
% is important to choose carefully the limits of the slope calculation (beg, fin)

%Choose the Q max
MaxQ = 10;
%It will sweep a symetric Q with 1000 points
Q = linspace(-MaxQ, MaxQ, 10);

% A square matrix is needed
[A,eps, tau, IQ] = MultifractalV2( Matriz, Q );
alpha = 0;
tau1 = 0;
figure(5);
%Limits for calculating the slope of TAU and A
beg = 1 ;
fin = 2;
hold on
 for i = 1:length(Q)


%Calculate the slope of tau, to que get Tau and A to get alpha
 p = polyfit(log(eps(beg:fin)'), log(tau(beg:fin,i)), 1);
tau1(i) = p(1);

p2 = polyfit(log(eps(beg:fin)'), log(A(beg:fin,i)), 1);
alpha(i) = p2(1);

%Uncomment to check if the chosen limits are correct

     figure(10)
    hold on
plot(log(eps(beg:fin)), log(tau(beg:fin,i)), '.-', 'MarkerSize', 30)


figure(11)
 hold on
  plot(log(eps(beg:fin)), log(A(beg:fin,i)), '.-', 'MarkerSize', 30)
 pause
 
% 


end

%Calculate Dq and Falpha
D = tau1./(Q-1);
F = alpha.*Q - tau1;
%Variables to export in origin for example
Save = [alpha', F'];
SaveQ = [Q', D'];

%Show Results
figure(25)
plot(Q, alpha)
hold on
figure(568)
hold on
plot(Q, tau1)
figure(56)
hold on
plot(alpha, F, 'o-')

%EE571_Introduction_to_Convex_Optimization_Homework 7

%Problem 3

%Given:

%The problem instance

%Part A : Implementing the 'illum_data' file to plot the lamp/patch geometry

% illum_data: generates the input data for the illumination problem 
%             (a matrix A whose rows are a_k^T)

figure(1)
L = [linspace(0,1,10); 
     1.9 1.8 1.0 1.1 1.9 1.8 1.9 1.7 1.5 1.5 ];  %lamp power 

m = size(L,2);    % number of lamps

% begin and endpoints of patches 
V = [linspace(0,1,21);
    .4* [0.0 0.1 0.15 0.2 0.1 0.2 0.3 0.0 0.0 0.0 ,  ...
      0.1 0.2 0.2 0.0 0.1 0.05 0.1 0.1 0.0 0.2 0.1]];
  
d=plot(L(1,:), L(2,:), 'bo', V(1,:), V(2,:), '-r');
title('The Problem Instance');

n = size(V,2)-1;  % number of patches

% construct A
dV = V(:,2:n+1)-V(:,1:n);    % tangent to patches
VI = V(:,1:n)+.5*dV;         % midpoint of patches
A = zeros(n,m);
for i=1:n
  for j=1:m
    dVI = L(:,j)-VI(:,i);  
    dVperp = null(dV(:,i)');  % upward pointing normal 
    if dVperp(2)<0
      dVperp = -dVperp;
    end
    A(i,j) = max(0,dVI'*dVperp/(norm(dVI)*norm(dVperp)))./norm(dVI)^2;
  end
end

%% Problem 3.1 Equal lamp powers

figure(2);
hold on;

%Choosing a value for gamma between 0 and 1
for gamma = 0:0.01:1 
    
    %Lamp power
    p01 = repmat(gamma,m,1);
    
    %Objective value for each gamma
    f01 = max(abs(log(A*p01)));
    
    disp(f01);
    
    grid on;
    plot(gamma,f01,'--bo');
    xlabel('Lamp power : p = gamma(0 <= gamma <= 1)');
    ylabel('Illumination produced by varying gamma : f0');
    title('Problem 3.1 Equal Lamp Powers');
    legend('Objective value for each gamma');
    
end

hold off;

%% Problem 3.2  Least Squares with saturation

    b02 = ones(n,1);
    
    %The optimal point for the Saturated Least Squares problem
    p02 = A \ b02;    
     
 %Comparing the values of the optimal solutions and setting them to zero or one    
 for i = 1: length(p02)   
    
     %If the value of the least square solution is negative, set it to zero 
     if(p02(i) < 0)
        p02(i) = 0;   
     
     %If the value of the least square solution is greater than one, set it to one   
     else if(p02(i) > 1)
             p02(i) = 1;
        end
    end
 end
 
%The optimal value for the Saturated Least Squares problem 
f02 = max(abs(log(A*p02)))        

%% Problem 3.3 Regularized Least Squares 

 b31 = ones(n,1);
 
 b32 = ones(m,1);

 %Varying the value of 'rho' such that the values in 'A' are in the range
 %of [0,1]
for rho = 0.1:0.01:10
   
   p31 = [A; sqrt(rho)*eye(m)]\[b31;sqrt(rho)*0.5*b32];
   
   t31 = zeros(m,1);
   t32 = ones(m,1);
    
   if(t31 <= p31 <= t32)
     
       break;
       
   end
end

%The optimal point for the Regularised Least Squares problem
p03 = p31;

%The optimal value for the Regularised Least Squares problem 
f03 = max(abs(log(A*p03)));

%% Problem 3.4 Chebyshev approximation

b04 = ones(n,1);
t41 = zeros(m,1);
t42 = ones(m,1);

cvx_begin
    
    variable p41(10)
        
        minimize (norm((A*p41) - b04,inf))
        subject to
                t41 <= p41 <= t42
                
cvx_end

%The optimal point for the Chebyshev Approximation problem
p4 = p41

%The optimal value for the Chebyshev Approximation problem
f04 = max(abs(log(A*p4)))

%% Problem 3.5 Exact solution

t51 = zeros(m,1);
t52 = ones(m,1);

cvx_begin
    
    variable p51(10)
        
        minimize (max([(A*p51);inv_pos(A*p51)]))
            subject to
                    t51 <= p51 <= t52
                    
cvx_end    

%The optimal point for the Chebyshev Approximation problem
p5 = p51

%The optimal value for the Chebyshev Approximation problem
f05 = max(abs(log(A*p5)))
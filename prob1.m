%EE571_Introduction_to_Convex_Optimization_Homework 7

%Problem 1

%Given:

%Resources consumed
A = [1 2 0 1;0 0 3 1;0 3 1 1;2 1 2 5;1 0 3 2];

%Maximum allowable limit for resource consumption 
cmax = repmat(100,5,1);

%Basic price
p = [3;2;7;6];

%Discount price per quantity
pdisc = [2;1;4;2];

%Quantity discount level
q = [4;10;5;10];

%% Method 1 : Formulating the given problem as a Linear Program

echo on

cvx_begin

    variable u(4)
    variable x(4)

    maximize (sum(1.'*u))
    subject to 
            x >= 0;
            (A * x <= cmax);
            (p.* x) >= u;
            ((p.*q)+(pdisc.*(x - q))) >= u;
cvx_end

%Display the activity levels
x

%Total revenue
r = min((p.*x),(p.*q)+(pdisc.*(x - q)))
total = sum(r)

%Revenue associated with each level
avg_price = r./x

echo off

%% Method 2 : Formulating the given problem as a Convex Problem

echo on

cvx_begin
    variable x(4)
    
    maximize (sum(min(p.*x,(p.*q + pdisc.*(x - q)))))
    subject to
            (A*x) <= cmax
            x >= 0
cvx_end

%Display the activity levels
x

%Total revenue
r = min((p.*x),(p.*q)+(pdisc.*(x - q)))
total = sum(r)

%Revenue associated with each level
avg_price = r./x
            
echo off

    
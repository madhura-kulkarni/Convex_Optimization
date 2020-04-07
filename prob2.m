%EE571_Introduction_to_Convex_Optimization_Homework 7

%Problem 2

%%(a) norm([x + 2y, x - y]) = = 0

%The equivalent reformulation is as follows:  
        x == 0;
        y == 0;
        
%% (b) square(square(x + y)) <= x - y

%The equivalent reformulation is as follows:  
        (x + y)^4 <= x - y;
        
%% (c) 1/x + 1/y <= 1; x >= 0; y >= 0

%The equivalent reformulation is as follows:  
            inv_pos(x) + inv_pos(y) <= 1;
           
%% (d) norm([max(x,1), max(y,2)]) <= 3*x + y

%The equivalent reformulation is as follows:  
            norm([u;v]) <= 3*x + y;
            max(x,1) <= u;
            max(y,2) <= v;
            
%% (e) x*y >= 1; x>= 0; y>= 0

%The equivalent reformulation is as follows:  
            x >= inv_pos(y);
            x >= 0;
            y >= 0;
            
%% (f) (x + y)^2/sqrt(y) <= x - y + 5

%The equivalent reformulation is as follows:  
    
            quad_over_lin(x + y,sqrt(y)) <= x - y + 5;
            
%% (g) x^3 + y^3 <= 1; x >= 0; y >= 0

%The equivalent reformulation is as follows:  

            pow_pos(x,3) + pow_pos(y,3) <= 1;
            
%% (h) x+z <= 1 + sqrt(x*y - z^2); x>= 0; y>= 0

%The equivalent reformulation is as follows:  

            x + z <= 1 + geo_mean([x - quad_over_lin(z,y),y]);
            
%% A convex problem solving these constraints:

echo on;

      cvx_begin
            variables x y u v z
            x == 0;
            y == 0;
            (x + y)^4 <= x - y;
            inv_pos(x) + inv_pos(y) <= 1;
            norm([u;v]) <= 3*x + y;
            max(x,1) <= u;
            max(y,2) <= v;
            x >= inv_pos(y);
            x >= 0;
            y >= 0;
            quad_over_lin(x + y,sqrt(y)) <= x - y + 5;
            pow_pos(x,3) + pow_pos(y,3) <= 1;
            x + z <= 1 + geo_mean([x-quad_over_lin(z,y),y]);
       cvx_end
       
echo off
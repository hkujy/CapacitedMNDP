% reference:
% The linear convergence of a successive linear programming algorithm
% Nonlinear optimization by successive linear programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the main function code
clc;clear;

%%% Step 1: solve the base case SUE model 
d01 = 100; 
d02 = 100;
A0 = [];
b0=[];
Aeq0 = [1,1,0,0,0;0,0,1,1,1];
beq0 = [d01;d02];
lb0 = [0;0;0;0;0];
ub0 = [inf;inf;inf;inf;inf];
m = 1;
l = 1;
for i=1:1:2
    x0(i) = d01/2;
end
for i=3:1:5
    x0(i) = d02/3;
end
x0 = fmincon(@(x)SUE_fun0(x), x0, A0, b0, Aeq0, beq0, lb0, ub0);
%%%

%%% set the parameters
% remark: Q max includes the initial demand d0
Q1 = 500;
Q2 = 500;
step_size = 30;
step_size_update_ratio = 0.3;
radius = 30;
%%% 

%%% set the capacity and test convergence 
for b52 = 100:100:100
    for b73 = 100:100:100
        fmax1 = Q1;
        fmax2 = Q2;
        d1(1) = 0;
        d1(2) = 100;
        d2(1) = 0;
        d2(2) = 100;
        n=2;
        iter = 0;
        % the following set the initial x1 values 
        % afterwards, x1 is set based on the x value from previous iterations
        x1(1) = d1(n) * 0.25;
        x1(2) = d1(n) * 0.25;
        x1(3) = d1(n) * 0.25;
        x1(4) = d1(n) * 0.25;
        x1(5) = d2(n) * 0.2;
        x1(6) = d2(n) * 0.2;
        x1(7) = d2(n) * 0.2;
        x1(8) = d2(n) * 0.2;
        x1(9) = d2(n) * 0.2;
        while ((abs(d1(n)-d1(n-1))/d1(n-1)>0.001||abs(d2(n)-d2(n-1))/d2(n-1)>0.001)&&iter<100)
            iter = iter+1;
            A1 = [];
            b1 = [];
            Aeq1 = [1,1,1,1,0,0,0,0,0;0,0,0,0,1,1,1,1,1];
            beq1 = [d1(n);d2(n)];
            lb1 = [0;0;0;0;0;0;0;0;0];
            ub1 = [inf;inf;inf;inf;inf;inf;inf;inf;inf];
            [x,fval,exitfalg,output,lambda(l),grad,hessian] = fmincon(@(x)SUE_fun(x,x0), x1, A1, b1, Aeq1, beq1, lb1, ub1);
            y = x + [x0(1),x0(2),0,0,x0(3),x0(4),x0(5),0,0];
            % del: path-link matrix
            del = [
                1,0,0,1,1,0,0,0,0,0,0,0;
                0,0,1,0,0,0,0,0,0,0,0,1;
                0,1,0,1,1,1,0,0,0,0,0,0;
                0,0,1,0,1,0,0,0,0,0,1,0;
                0,0,0,0,1,0,0,1,0,0,0,0;
                0,0,0,0,0,0,0,0,0,1,0,1;
                0,0,0,1,1,0,1,0,0,0,0,0;
                0,0,0,1,1,1,0,0,1,0,0,0;
                0,0,0,0,1,0,0,0,0,1,1,0];
            v = y*del;
            o1 = y(1)+y(2)+y(3)+y(4);
            o2 = y(5)+y(6)+y(7)+y(8)+y(9);
            f1(1) = 1;
            f1(2) = 1;
            od_path = [1 1 1 1 0 0 0 0 0];
            for i=1:1:12
                f21(i) = y.* od_path* del(:,i)/(d01+d1(n));
            end
            od_path = [0 0 0 0 1 1 1 1 1];
            for i=1:1:12
                f22(i) = y.* od_path*del(:,i)/(d02+d2(n)); 
            end
            
            A2 = [
                f21(6),f22(6);     
                f21(11),f22(11);
                f1(1),0;
                0,f1(2)
                ];
            s6 = b52-v(6)+f21(6)*d1(n)+f22(6)*d2(n);
            s11 = b73-v(11)+f21(11)*d1(n)+f22(11)*d2(n);
            % update the boundary
            this_iter_max1=min(d1(n)+step_size,fmax1);
            this_iter_max2=min(d2(n)+step_size,fmax2);
            this_iter_min1=max(d1(n)-step_size,0);
            this_iter_min2=max(d2(n)-step_size,0);

            b2 = [
                s6;              
                s11;              
                this_iter_max1+d01-(y(1)+y(2)+y(3)+y(4))+f1(1)*d1(n);
                this_iter_max2+d02-(y(5)+y(6)+y(7)+y(8)+y(9))+f1(2)*d2(n)
                ];
            Aeq2 = [];
            beq2 = [];
            lb2 = [this_iter_min1;this_iter_min2];
            ub2 = [this_iter_max1;this_iter_max2];
            fun = [-1,-1];
            n = n+1;
            options = optimoptions(@linprog,'Display','iter');
            [d,Lp_fval,Lp_exitflag,Lp_output,Lp_lambda] = linprog(fun,A2,b2,Aeq2,beq2,lb2,ub2,options);
            d1(n) = d(1,1);
            d2(n) = d(2,1);
            fprintf('Gap = %10.4f,',max(abs(d1(n)-d1(n-1))/d1(n-1),abs(d2(n)-d2(n-1))/d2(n-1)))
            fprintf('d1 = %10.4f, d2 = %10.4f\n',d1(n),d2(n))
            fprintf('current: radius gap = %10.4f, step_size = %10.4f\n',abs(d1(n)+d2(n)-d1(n-1)-d2(n-1)), step_size)
            if (abs(d1(n)+d2(n)-d1(n-1)-d2(n-1))<0.67*radius)
                step_size = step_size*step_size_update_ratio;
                radius =  radius *step_size_update_ratio;
                fprintf('update step_size: step_size = %10.4f, radius = %10.4f\n',step_size,radius)
            end 
            % the following set the initial flow ratio for the next iteration
            x1(1)=d1(n)*(x(1)/d1(n-1));
            x1(2)=d1(n)*(x(2)/d1(n-1));
            x1(3)=d1(n)*(x(3)/d1(n-1));
            x1(4)=d1(n)*(x(4)/d1(n-1));
            x1(5)=d2(n)*(x(5)/d2(n-1));
            x1(6)=d2(n)*(x(6)/d2(n-1));
            x1(7)=d2(n)*(x(7)/d2(n-1));
            x1(8)=d2(n)*(x(8)/d2(n-1));
            x1(9)=d2(n)*(x(9)/d2(n-1));           
            y = x + [x0(1),x0(2),0,0,x0(3),x0(4),x0(5),0,0]; 
            v = y*del;
            fprintf('Check Cap: itera = %d,Link 6 Flow = %10.4f,  Link 11 Flow = %10.4f\n',iter,v(6), v(11))
        end
        % the following compute and double check the final solution value
        x1 = x;
        [x,fval,exitfalg,output,lambda(l),grad,hessian] = fmincon(@(x)SUE_fun(x,x0), x1, A1, b1, Aeq1, beq1, lb1, ub1);
        y = x + [x0(1),x0(2),0,0,x0(3),x0(4),x0(5),0,0];
        v = y*del;
        fprintf('Final Check Cap: Link 6 Flow = %10.4f,  Link 11 Flow = %10.4f\n',v(6), v(11))

        p{m} = y;
        
        m = m+1;
        
        o(l) = d1(n) + d2(n);
        
        t0=[20 5 5 3 21 2 12 25 3 4 3 40]; 
        
        t(1) = t0(1);
        t(2) = t0(2);
        t(3) = t0(3)+ 0.15*(v(3)/500)*(v(3)/500)*(v(3)/500)*(v(3)/500);
        t(4) = t0(4)+ v(4)/500;
        t(5) = t0(5)+ v(5)/500;
        t(6) = t0(6);
        t(7) = t0(7);
        t(8) = t0(8);
        t(9) = t0(9);
        t(10) = t0(10)+0.15*(v(10)/500)*(v(10)/500)*(v(10)/500)*(v(10)/500);
        t(11) = t0(11);
        t(12) = t0(12)+0.15*(v(12)/1000)*(v(12)/1000)*(v(12)/1000)*(v(12)/1000);
        
        tp1(l) = t(1) + t(4) + t(5);
        tp2(l) = t(3) + t(12);
        tp3(l) = t(2) + t(6) + t(4) + t(5);
        tp4(l) = t(3) + t(11) + t(5);
        tp5(l) = t(8) + t(5);
        tp6(l) = t(10) + t(12);
        tp7(l) = t(7) + t(4) + t(5);
        tp8(l) = t(9) + t(6) + t(4) + t(5);
        tp9(l) = t(10) + t(11) + t(5);
        
        tp0(1) = t0(1) + t0(4) + t0(5);
        tp0(2) = t0(3) + t0(12);
        tp0(3) = t0(2) + t0(6) + t0(4) + t0(5);
        tp0(4) = t0(3) + t0(11) + t0(5);
        tp0(5) = t0(8) + t0(5);
        tp0(6) = t0(10) + t0(12);
        tp0(7) = t0(7) + t0(4) + t0(5);
        tp0(8) = t0(9) + t0(6) + t0(4) + t0(5);
        tp0(9) = t0(10) + t0(11) + t0(5);
        
        tc(l) = tp1(l)*y(1) + tp2(l)*y(2) + tp3(l)*y(3) + tp4(l)*y(4) + tp5(l)*y(5) + tp6(l)*y(6) + tp7(l)*y(7) + tp8(l)*y(8) + tp9(l)*y(9);
        
        l = l+1;
    end
    
end
save final

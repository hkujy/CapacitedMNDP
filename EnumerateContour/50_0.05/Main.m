% reference:
% The linear convergence of a successive linear programming algorithm
% Nonlinear optimization by successive linear programming
clc;clear;

d01 = 100;
d02 = 100;
maxIter = 200;
epsGap = 0.001;

A0 = [];
b0=[];
Aeq0 = [1,1,0,0,0;0,0,1,1,1];
beq0 = [d01;d02];
lb0  = [0;0;0;0;0];
ub0  = [d01;d01;d02;d02;d02];

m = 1;
l = 1;
for i=1:1:2
    x0(i) = d01/2;
end

for i=3:1:5
    x0(i) = d02/3;
end

theta = 0.05;
x0 = fmincon(@(x)SUE_fun0(x,theta), x0, A0, b0, Aeq0, beq0, lb0, ub0);
%save x0;
%load x0;
% remark: Q max includes the intial demand d0
Q1 = 1000;
Q2 = 1000;

b52_iniCap = 0;
b52_incre = 50;
b52_maxCap = 400;
b52_num = (b52_maxCap-b52_iniCap)/b52_incre + 1;

b73_iniCap = 0;
b73_incre = 50;
b73_maxCap = 400;
b73_num = (b73_maxCap-b73_iniCap)/b73_incre + 1;

CaseIndex = 0;
fileID_1 = fopen('DemandData.txt','w');
fileID_2 = fopen('PathData.txt','w');
fileID_3 = fopen('LinkData.txt','w');

fprintf(fileID_1, 'CaseId,b52,b73,d1,d2,sum\n');
fprintf(fileID_2, 'CaseId,b52,b73,pid,flow,cost\n');
fprintf(fileID_3, 'CaseId,b52,b73,lid,flow,cost\n');

for b52 = b52_iniCap:b52_incre:b52_maxCap  %bike and ride
    for b73 = b73_iniCap:b73_incre:b73_maxCap  %park and ride
        CaseIndex = CaseIndex + 1;
        fmax1 = Q1;
        fmax2 = Q2;
        d1(1) = 0;
        d1(2) = 30;
        d2(1) = 0;
        d2(2) = 30;
        n=2;
        iter = 0;
        minCap = 30;
        if (b52>0&&b73>0)
            minCap = min(b52,b73);
	        del52 = 1;
            del73 = 1;
            ca52 = max(Q1,Q2);
            ca73 = max(Q1,Q2);
        end
        if (b52>0&&b73==0)
            minCap = b52;
            del52 = 1;
            ca52 = max(Q1,Q2);
            del73 = 0;
            ca73 = 0;
        end
        if (b52==0&&b73>0)
            minCap = b73;
            del52 = 0;
	        del73 = 1;
	        ca52 = 0;
	        ca73 = max(Q1,Q2);
        end
        step_size(1) = min(30,minCap);
        step_size(2) = min(30,minCap);
        radius(1) = 30;
        radius(2) = 30;
        step_size_update_ratio = 0.6;
        
        % the following set the initial x1 values 
        % afterwards, x1 is set based on the x value from previous iterations
        x1(1) = d1(n);
        x1(2) = 0;
        x1(3) = 0;
        x1(4) = 0;
        x1(5) = d2(n);
        x1(6) = 0;
        x1(7) = 0;
        x1(8) = 0;
        x1(9) = 0;
 
        % remark: exist mulitple solutions
        while ((abs(d1(n)-d1(n-1))/d1(n-1)>epsGap||abs(d2(n)-d2(n-1))/d2(n-1)>epsGap)...
                &&(b52>0||b73>0)&&iter<maxIter)

            iter = iter+1;
            A1 = [];
            b1 = [];
            Aeq1 = [1,1,del52,del73,0,0,0,0,0;
                    0,0,0,0,1,1,1,del52,del73];
            beq1 = [d1(n);d2(n)];
            lb1 = [0;0;0;0;0;0;0;0;0];
            ub1 = [Q1;Q1;ca52;ca73;Q2;Q2;Q2;ca52;ca73];
           
            [x,fval,exitfalg,output,lambda(l),grad,hessian] = fmincon(@(x)SUE_fun(x,x0,b52,b73,theta), x1, A1, b1, Aeq1, beq1, lb1, ub1);
            y = x + [x0(1),x0(2),0,0,x0(3),x0(4),x0(5),0,0];
            
            % del: path-link matrix
            del = [
                1,0,0,1,1,0,0,0,0,0,0,0;
                0,0,1,0,0,0,0,0,0,0,0,1;
                0,del52,0,del52,del52,del52,0,0,0,0,0,0;
                0,0,del73,0,del73,0,0,0,0,0,del73,0;
                0,0,0,0,1,0,0,1,0,0,0,0;
                0,0,0,0,0,0,0,0,0,1,0,1;
                0,0,0,1,1,0,1,0,0,0,0,0;
                0,0,0,del52,del52,del52,0,0,del52,0,0,0;
                0,0,0,0,del73,0,0,0,0,del73,del73,0];
            v = y*del;

            o1 = y(1)+y(2)+y(3)+y(4);
            o2 = y(5)+y(6)+y(7)+y(8)+y(9);
            
            f1(1) = 1;
            f1(2) = 1;
           % od_path incidence for OD pair 1
            od_path = [1 1 1 1 0 0 0 0 0];
            for i=1:1:12
                f21(i) = y.* od_path* del(:,i)/(d01+d1(n)); %#ok<*SAGROW
            end
           % od_path incidence for OD pair 2
            od_path = [0 0 0 0 1 1 1 1 1];
            for i=1:1:12
                f22(i) = y.* od_path*del(:,i)/(d02+d2(n)); %#ok<*SAGROW>        
            end
            
            A2 = [
                f21(6),f22(6);     
                f21(11),f22(11);
                f1(1),0;
                0,f1(2)
                ];
            s6 = b52-v(6)+f21(6)*d1(n)+f22(6)*d2(n);
            s11 = b73-v(11)+f21(11)*d1(n)+f22(11)*d2(n);
            % the following code update the upper and lower bounds       
            this_iter_max1=min(d1(n)+step_size(1),fmax1);
            this_iter_max2=min(d2(n)+step_size(2),fmax2);
            this_iter_min1=max(d1(n)-5*step_size(1),5);
            this_iter_min2=max(d2(n)-5*step_size(2),5);
             
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
            fprintf('b52 = %4.1f, b73 = %4.1f, Iter = %3.0f, Gap = %6.4f,',b52,b73,iter,max(abs(d1(n)-d1(n-1))/d1(n-1),abs(d2(n)-d2(n-1))/d2(n-1)))
            fprintf('d1 = %6.2f, d2 = %6.2f\n',d1(n),d2(n))
            fprintf('current: radius gap = %6.2f, step_size(1) = %6.2f, step_size(2) = %6.2f\n',...
                    abs(d1(n)+d2(n)-d1(n-1)-d2(n-1)), step_size(1),step_size(2))

            if (abs(d1(n)-d1(n-1))<radius(1))
                step_size(1) = step_size(1)*step_size_update_ratio;
                radius(1) =  radius(1)*step_size_update_ratio;
                fprintf('update step_size(1): step_size = %10.4f, radius = %10.4f\n',step_size(1),radius(1))
            end
            if (abs(d2(n)-d2(n-1))< radius(2))
                step_size(2) = step_size(2)*step_size_update_ratio;
                radius(2) =  radius(2)*step_size_update_ratio;
                fprintf('update step_size(2): step_size = %10.4f, radius = %10.4f\n',step_size(2),radius(2))
            end


            % if (abs(d1(n)+d2(n)-d1(n-1)-d2(n-1))<0.67*radius)
            %     step_size = step_size*step_size_update_ratio;
            %     radius =  radius *step_size_update_ratio;
            %     fprintf('update step_size: step_size = %10.4f, radius = %10.4f\n',step_size,radius)
            % end 
            
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
        
        save(strcat('b52_',int2str(b52),'_b73_',int2str(b73)));
        
        if (b52==0&&b73==0)
            y = [x0(1),x0(2),0,0,x0(3),x0(4),x0(5),0,0];
            del = [
                1,0,0,1,1,0,0,0,0,0,0,0;
                0,0,1,0,0,0,0,0,0,0,0,1;
                0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,1,0,0,1,0,0,0,0;
                0,0,0,0,0,0,0,0,0,1,0,1;
                0,0,0,1,1,0,1,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0];
            v = y * del;
        end
        
        p{l} = y;
        
        m = m+1;
                
        o(l) = sum(y);
        
        % the following compute and double check the final solution value\
        if (b52>0||b73>0)
            x1 = x;
            [x,fval,exitfalg,output,lambda(l),grad,hessian] = fmincon(@(x)SUE_fun(x,x0,b52,b73,theta), x1, A1, b1, Aeq1, beq1, lb1, ub1);
            y = x + [x0(1),x0(2),0,0,x0(3),x0(4),x0(5),0,0];
            v = y*del;
            fprintf('Final Check Cap: Link 6 Flow = %10.4f,  Link 11 Flow = %10.4f\n',v(6), v(11))
        end
        
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
        t(10) = t0(10) + 0.15*(v(10)/500)*(v(10)/500)*(v(10)/500)*(v(10)/500);
        t(11) = t0(11);
        t(12) = t0(12) + 0.15*(v(12)/1000)*(v(12)/1000)*(v(12)/1000)*(v(12)/1000);
        

        tp1(l) = t(1) + t(4) + t(5);
        tp2(l) = t(3) + t(12);
        tp3(l) = t(2) + t(6) + t(4) + t(5);
        tp4(l) = t(3) + t(11) + t(5);
        tp5(l) = t(8) + t(5);
        tp6(l) = t(10) + t(12);
        tp7(l) = t(7) + t(4) + t(5);
        tp8(l) = t(9) + t(6) + t(4) + t(5);
        tp9(l) = t(10) + t(11) + t(5);

        tp_cost = [tp1(l),tp2(l),tp3(l),tp4(l),tp5(l),tp6(l),tp7(l),tp8(l),tp9(l)];
        
        
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
        % print demand cases
        fprintf(fileID_1, '%2.0f,%4.0f,%4.0f,%8.4f,%8.4f,%8.4f\n',...
                            CaseIndex,b52,b73,d1(n),d2(n),d1(n)+d2(n));
        for i = 1:9  % print path flow and cost
            fprintf(fileID_2, '%2.0f,%4.0f,%4.0f,%2.0f,%8.4f,%8.4f\n',...
                                CaseIndex,b52,b73,i,y(i),tp_cost(i));
        end
        for i = 1:12    % print link flow cost
            fprintf(fileID_3, '%2.0f,%4.0f,%4.0f,%2.0f,%8.4f,%8.4f\n',...
                                CaseIndex,b52,b73,i,v(i),t(i));
        end 
        

    end
    
end
save final

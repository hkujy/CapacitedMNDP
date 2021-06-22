function z = SUE_fun0(x,theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t0=[20 5 3 21 12 25 4 40];
syms a b c d e f g h;
del = [1,0,1,1,0,0,0,0;
       0,1,0,0,0,0,0,1;
       0,0,0,1,0,1,0,0;
       0,0,0,0,0,0,1,1;
       0,0,1,1,1,0,0,0];
   
v = x * del;

z1(1) = int(t0(1),a,0,v(1));
z1(2) = int(t0(2)+0.15*(b/500)*(b/500)*(b/500)*(b/500),b,0,v(2));
z1(3) = int(t0(3)+ c/500,c,0,v(3));
z1(4) = int(t0(4)+ d/500,d,0,v(4));
z1(5) = int(t0(5),e,0,v(5));
z1(6) = int(t0(6),f,0,v(6));
z1(7) = int(t0(7)+0.15*(g/500)*(g/500)*(g/500)*(g/500),g,0,v(7));
z1(8) = int(t0(8)+0.15*(h/1000)*(h/1000)*(h/1000)*(h/1000),h,0,v(8));

for l=1:1:5
    z2(l) = x(l)*(log(x(l))-1)/theta; %#ok<*AGROW>
end

o(1) = x(1)+x(2);
o(2) = x(3)+x(4)+x(5);
m(1) = x(1);
m(2) = x(3)+x(5);
c(1) = x(2);
c(2) = x(4);

for i=1:1:2    
    z41(i) = (m(i)*(log(m(i))-1))/theta;
    z42(i) = (c(i)*(log(c(i))-1))/theta;
end

z=double(sum(z1)+sum(z2)+sum(z41)+sum(z42));

end

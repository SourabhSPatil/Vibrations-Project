%% single mode
clc
clear
close
rhol = 1;
El = 1;
Il = 1;
Al = 1;
l=1;
Cb = [0,0.2,0.4,0.6,0.8];
Ch = 0;
mode_shape_coeff = ones(3,3,5);
wn = ones(3,5);
r = ones(3,5);
for i = 1:5
syms x w C1 C2
Ex = 1+(x/l);
rhox = 1+(x/l)+(x/l)^2;
Ix = (1-(Cb(i)*x/l))*(1-(Ch*x/l)^3);
Ax = (1-(Cb(i)*x/l))*(1-(Ch*x/l));
beta = (w^2 * rhol*Al/(El*Il))^(0.25);
function_coeff = -C1*(sin(pi*x/l));
double_diff = diff(function_coeff,2);
bracket = Ex*Ix*double_diff;
fourth_diff = diff(bracket,2);
first_term = fourth_diff*function_coeff/C1;
A = int(first_term,x,0,l);
second_term = beta^4 *rhox*Ax * function_coeff*function_coeff/C1;
B = int(second_term,x,0,l);
Z = w^2 *(A/B);
wn1(i) = sqrt(double(Z));
X = linspace(0,l,100);
Y = -(sin(pi*X/l));
plot(X,Y)
hold on
end
%% Three modes
clc
clear
close
rhol = 1;
El = 1;
Il = 1;
Al = 1;
l=1;
Cb = [0,0.2,0.4,0.6,0.8];
Ch = 0;
wn1(5) = 0;
for i = 1:5
for j =1:3
syms x w C1 C2 C3
Ex = 1+(x/l);
rhox = 1+(x/l)+(x/l)^2;
Ix = (1-(Cb(i)*x/l))*(1-(Ch*x/l)^3);
Ax = (1-(Cb(i)*x/l))*(1-(Ch*x/l));
beta = (w * rhol*Al/(El*Il))^(0.25);
function_coeff = [C1*(sin(pi*x/l)*(x-x^2)) C2*(sin(2*pi*x/l)*(x-x^2)) C3*(sin(3*pi*x/l)*(x-x^2))];
function_mat = [(sin(pi*x/l)*(x-x^2)) (sin(2*pi*x/l)*(x-x^2)) (sin(3*pi*x/l)*(x-x^2))];
constant = [C1 C2 C3];
for k = 1:3
double_diff = diff(function_coeff(k),2);
bracket = Ex*Ix*double_diff;
fourth_diff = diff(bracket,2);
first_term = fourth_diff*function_mat(j);
integ1 = int(first_term,x,0,l);
A(j,k) = double(integ1/constant(k));
end
for k = 1:3
second_term = -beta^4 *rhox*Ax * function_coeff(k)*function_mat(j);
integ2 = int(second_term,x,0,l);
B(j,k) = double(integ2/(constant(k)*w));
end
end
B_new = B*w;
final_mat = A + B*w;
determinant_exp = det(final_mat,'Algorithm','minor-expansion');
c = sym2poly(determinant_exp);
r(:,i) = sort(roots(c));
wn(:,i) = sqrt(r(:,i));
B1 = double(subs(B_new,w,r(1,i)));
B2 = double(subs(B_new,w,r(2,i)));
B3 = double(subs(B_new,w,r(3,i)));
final_mat1 = A+B1;
final_mat2 = A+B2;
final_mat3 = A+B3;
v1 = [0-final_mat1(2,1);0-final_mat1(3,1)];
v2 = [0-final_mat2(2,1);0-final_mat2(3,1)];
v3 = [0-final_mat3(2,1);0-final_mat3(3,1)];
mode_shape_coeff_single(1,i) = 1;
mode_shape_coeff_double(1,i) = 1;
mode_shape_coeff_triple(1,i) = 1;
mode_shape_coeff_single(2:3,i) = final_mat1(2:3,2:3)\v1;
mode_shape_coeff_double(2:3,i) = final_mat2(2:3,2:3)\v2;
mode_shape_coeff_triple(2:3,i) = final_mat3(2:3,2:3)\v3;
end
%% plots
X = linspace(0,l,100);
figure(1)
for i = 1:5
Y = mode_shape_coeff_single(1,i)*(sin(pi*X/l).*(X-X.^2))+ mode_shape_coeff_single(2,i)*(sin(2*pi*X/l).*(X-X.^2))+ mode_shape_coeff_single(3,i).*(sin(3*pi*X/l).*(X-X.^2));
plot(X,Y)
hold on
end
title("First Mode")
xlabel("x/l")
ylabel("response")
Legend = cell(5,1);
Legend{1} = num2str(Cb(1));
Legend{2} = num2str(Cb(2));
Legend{3} = num2str(Cb(3));
Legend{4} = num2str(Cb(4));
Legend{5} = num2str(Cb(5));
legend(Legend);
figure(2)
for i = 1:4
Y = mode_shape_coeff_double(1,i)*(sin(pi*X/l).*(X-X.^2))+ mode_shape_coeff_double(2,i)*(sin(2*pi*X/l).*(X-X.^2))+ mode_shape_coeff_double(3,i).*(sin(3*pi*X/l).*(X-X.^2));
plot(X,Y)
hold on
end
title("Second Mode")
xlabel("x/l")
ylabel("response")
Legend = cell(4,1);
Legend{1} = num2str(Cb(1));
Legend{2} = num2str(Cb(2));
Legend{3} = num2str(Cb(3));
Legend{4} = num2str(Cb(4));
legend(Legend);
figure(3)
for i = 1:5
Y = mode_shape_coeff_triple(1,i)*(sin(pi*X/l).*(X-X.^2))+ mode_shape_coeff_triple(2,i)*(sin(2*pi*X/l).*(X-X.^2))+ mode_shape_coeff_triple(3,i).*(sin(3*pi*X/l).*(X-X.^2));
plot(X,Y)
hold on
end
title("Third Mode")
xlabel("x/l")
ylabel("response")
Legend = cell(5,1);
Legend{1} = num2str(Cb(1));
Legend{2} = num2str(Cb(2));
Legend{3} = num2str(Cb(3));
Legend{4} = num2str(Cb(4));
Legend{5} = num2str(Cb(5));
legend(Legend);




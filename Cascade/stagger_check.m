clear 
close all
load cascade_fp_stagger5;
L_unsteady1_orig = L_unsteady;
h = 0.3;
Lift_Panel_Interp_4000steps;
L_unsteady1 = L_unsteady_fix;
foft = L_unsteady;
ffer3;
close figure 1
fer_Panel = fer2;
f_Panel = f;
figure(11112)
hp5 = semilogy(f_Panel*pi,abs(fer_Panel/vortstrength)./exp(-f_Panel*h*pi/1)*2*pi,'o');
hold on;


% Import relevant data
filename = 'ventres_fp_stagger5.txt';
A = importdata(filename);

for i = 1:length(A)
    k_ven(i) = A(i,1);
    cl_ven(i) = A(i,2);
   
end
cl_ven = cl_ven%.*exp(-k_ven*h);
figure(11112)
hp = semilogy(k_ven,cl_ven);
xlim([0 20])
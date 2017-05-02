% Final Project Matlab Plotting 
% Written By: Robert Masti
% This function will read in the data and serve as a template for the different cases
clc; close all; clear all;
num_ghost=3;
NEQ=4;

%% Mesh 1
x9 = load('Case_1_Mesh_1_Flux_1_order_2/xc-0.txt');
y9 = load('Case_1_Mesh_1_Flux_1_order_2/yc-0.txt');
ni9 = length(x9(1,:)); 
nj9 = length(x9(:,1));

mms9 = zeros(nj9, ni9, NEQ);
num9 = zeros(nj9, ni9, NEQ);
% conserved vars
rhou9 = load('Case_1_Mesh_1_Flux_1_order_2/rhou-79.txt');
rhov9 = load('Case_1_Mesh_1_Flux_1_order_2/rhov-79.txt');
rhoe9 = load('Case_1_Mesh_1_Flux_1_order_2/rhoe-79.txt');
% primvar mms
temp = load('Case_1_Mesh_1_Flux_1_order_2/rho_MMS-0.txt');
mms9(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_order_2/u_MMS-0.txt');
mms9(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_order_2/v_MMS-0.txt');
mms9(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_order_2/p_MMS-0.txt');
mms9(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_1_Flux_1_order_2/rho-79.txt');
num9(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_order_2/u-79.txt');
num9(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_order_2/v-79.txt');
num9(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_order_2/p-79.txt');
num9(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);


%% Mesh 2
x17 = load('Case_1_Mesh_2_Flux_1_order_2/xc-0.txt');
y17 = load('Case_1_Mesh_2_Flux_1_order_2/yc-0.txt');
ni17 = length(x17(1,:)); 
nj17 = length(x17(:,1));
mms17 = zeros(nj17, ni17, NEQ);
num17 = zeros(nj17, ni17, NEQ);
% conserved vars
rhou17 = load('Case_1_Mesh_2_Flux_1_order_2/rhou-407.txt');
rhov17 = load('Case_1_Mesh_2_Flux_1_order_2/rhov-407.txt');
rhoe17 = load('Case_1_Mesh_2_Flux_1_order_2/rhoe-407.txt');
% primvar mms
temp = load('Case_1_Mesh_2_Flux_1_order_2/rho_MMS-0.txt');
mms17(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_order_2/u_MMS-0.txt');
mms17(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_order_2/v_MMS-0.txt');
mms17(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_order_2/p_MMS-0.txt');
mms17(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_2_Flux_1_order_2/rho-407.txt');
num17(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_order_2/u-407.txt');
num17(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_order_2/v-407.txt');
num17(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_order_2/p-407.txt');
num17(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 3
x33 = load('Case_1_Mesh_3_Flux_1_order_2/xc-0.txt');
y33 = load('Case_1_Mesh_3_Flux_1_order_2/yc-0.txt');
ni33 = length(x33(1,:)); 
nj33 = length(x33(:,1));
mms33 = zeros(nj33, ni33, NEQ);
num33 = zeros(nj33, ni33, NEQ);
% conserved vars
rhou33 = load('Case_1_Mesh_3_Flux_1_order_2/rhou-2648.txt');
rhov33 = load('Case_1_Mesh_3_Flux_1_order_2/rhov-2648.txt');
rhoe33 = load('Case_1_Mesh_3_Flux_1_order_2/rhoe-2648.txt');
% primvar mms
temp = load('Case_1_Mesh_3_Flux_1_order_2/rho_MMS-0.txt');
mms33(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_order_2/u_MMS-0.txt');
mms33(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_order_2/v_MMS-0.txt');
mms33(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_order_2/p_MMS-0.txt');
mms33(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_3_Flux_1_order_2/rho-2648.txt');
num33(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_order_2/u-2648.txt');
num33(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_order_2/v-2648.txt');
num33(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_order_2/p-2648.txt');
num33(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 4
x65 = load('Case_1_Mesh_4_Flux_1_order_2/xc-0.txt');
y65 = load('Case_1_Mesh_4_Flux_1_order_2/yc-0.txt');
ni65 = length(x65(1,:)); 
nj65 = length(x65(:,1));
mms65 = zeros(nj65, ni65, NEQ);
num65 = zeros(nj65, ni65, NEQ);
% conserved vars
rhou65 = load('Case_1_Mesh_4_Flux_1_order_2/rhou-638.txt');
rhov65 = load('Case_1_Mesh_4_Flux_1_order_2/rhov-638.txt');
rhoe65 = load('Case_1_Mesh_4_Flux_1_order_2/rhoe-638.txt');
% primvar mms
temp = load('Case_1_Mesh_4_Flux_1_order_2/rho_MMS-0.txt');
mms65(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_order_2/u_MMS-0.txt');
mms65(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_order_2/v_MMS-0.txt');
mms65(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_order_2/p_MMS-0.txt');
mms65(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_4_Flux_1_order_2/rho-638.txt');
num65(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_order_2/u-638.txt');
num65(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_order_2/v-638.txt');
num65(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_order_2/p-638.txt');
num65(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 5
x129 = load('Case_1_Mesh_5_Flux_1_order_2/xc-0.txt');
y129 = load('Case_1_Mesh_5_Flux_1_order_2/yc-0.txt');
ni129 = length(x129(1,:)); 
nj129 = length(x129(:,1));
mms129 = zeros(nj129, ni129, NEQ);
num129 = zeros(nj129, ni129, NEQ);
% conserved vars
rhou129 = load('Case_1_Mesh_5_Flux_1_order_2/rhou-2750.txt');
rhov129 = load('Case_1_Mesh_5_Flux_1_order_2/rhov-2750.txt');
rhoe129 = load('Case_1_Mesh_5_Flux_1_order_2/rhoe-2750.txt');
% primvar mms
temp = load('Case_1_Mesh_5_Flux_1_order_2/rho_MMS-0.txt');
mms129(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_order_2/u_MMS-0.txt');
mms129(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_order_2/v_MMS-0.txt');
mms129(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_order_2/p_MMS-0.txt');
mms129(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_5_Flux_1_order_2/rho-2750.txt');
num129(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_order_2/u-2750.txt');
num129(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_order_2/v-2750.txt');
num129(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_order_2/p-2750.txt');
num129(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 6
x257 = load('Case_1_Mesh_6_Flux_1_order_2/xc-0.txt');
y257 = load('Case_1_Mesh_6_Flux_1_order_2/yc-0.txt');
ni257 = length(x257(1,:)); 
nj257 = length(x257(:,1));
mms257 = zeros(nj257, ni257, NEQ);
num257 = zeros(nj257, ni257, NEQ);
% conserved vars
rhou257 = load('Case_1_Mesh_6_Flux_1_order_2/rhou-6523.txt');
rhov257 = load('Case_1_Mesh_6_Flux_1_order_2/rhov-6523.txt');
rhoe257 = load('Case_1_Mesh_6_Flux_1_order_2/rhoe-6523.txt');
% primvar mms
temp = load('Case_1_Mesh_6_Flux_1_order_2/rho_MMS-0.txt');
mms257(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_order_2/u_MMS-0.txt');
mms257(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_order_2/v_MMS-0.txt');
mms257(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_order_2/p_MMS-0.txt');
mms257(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_6_Flux_1_order_2/rho-6523.txt');
num257(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_order_2/u-6523.txt');
num257(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_order_2/v-6523.txt');
num257(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_order_2/p-6523.txt');
num257(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Get the Error's

Error9 = num9-mms9;
Error17 = num17-mms17;
Error33 = num33-mms33;
Error65 = num65-mms65;
Error129 = num129-mms129;
Error257 = num257-mms257;

RelError9 = Error9./mms9;
RelError17 = Error17./mms17;
RelError33 = Error33./mms33;
RelError65 = Error65./mms65;
RelError129 = Error129./mms129;
RelError257 = Error257./mms257;

Linfnorm = zeros(6,NEQ);L2norm = zeros(6,NEQ);L1norm = zeros(6,NEQ);
for eq = 1:NEQ
    L2norm(1,eq) = sqrt(norm(Error9(:,:,eq),2)^2/(ni9*nj9));
    L2norm(2,eq) = sqrt(norm(Error17(:,:,eq),2)^2/(ni17*nj17));
    L2norm(3,eq) = sqrt(norm(Error33(:,:,eq),2)^2/(ni33*nj33));
    L2norm(4,eq) = sqrt(norm(Error65(:,:,eq),2)^2/(ni65*nj65));
    L2norm(5,eq) = sqrt(norm(Error129(:,:,eq),2)^2/(ni129*nj129));
    L2norm(6,eq) = sqrt(norm(Error257(:,:,eq),2)^2/(ni257*nj257));
        
    L1norm(1,eq) = (norm(Error9(:,:,eq),1)/(ni9*nj9));
    L1norm(2,eq) = (norm(Error17(:,:,eq),1)/(ni17*nj17));
    L1norm(3,eq) = (norm(Error33(:,:,eq),1)/(ni33*nj33));
    L1norm(4,eq) = (norm(Error65(:,:,eq),1)/(ni65*nj65));
    L1norm(5,eq) = (norm(Error129(:,:,eq),1)/(ni129*nj129));
    L1norm(6,eq) = (norm(Error257(:,:,eq),1)/(ni257*nj257));
    
    Linfnorm(1,eq) = norm(Error9(:,:,eq),inf);
    Linfnorm(1,eq) = norm(Error17(:,:,eq),inf);
    Linfnorm(1,eq) = norm(Error33(:,:,eq),inf);
    Linfnorm(1,eq) = norm(Error65(:,:,eq),inf);
    Linfnorm(1,eq) = norm(Error129(:,:,eq),inf);
    Linfnorm(1,eq) = norm(Error257(:,:,eq),inf);
    
end


% Plotting
h = [1,2,4,8,16];

OOA = zeros(length(h),NEQ,3);

for i = 1:length(h)
    for eq=1:NEQ
        OOA(i,eq,1) = abs(log((L2norm(i+1,eq)/L2norm(i,eq)))/log(2));         
        OOA(i,eq,2) = abs(log((L1norm(i+1,eq)/L1norm(i,eq)))/log(2));             
        OOA(i,eq,3) = abs(log((Linfnorm(i+1,eq)/Linfnorm(i,eq)))/log(2));       
    end
end

figure(1)
% L2 norm
semilogx(h,OOA(:,1,1),'r-^', 'MarkerFaceColor','r','MarkerSize',6)
hold on
semilogx(h,OOA(:,2,1),'b-^', 'MarkerFaceColor','b','MarkerSize',6)
semilogx(h,OOA(:,3,1),'g-^', 'MarkerFaceColor','g','MarkerSize',6)
semilogx(h,OOA(:,4,1),'k-^', 'MarkerFaceColor','k','MarkerSize',6)
% Linf norm
semilogx(h,OOA(:,1,2),'r--d', 'MarkerFaceColor','r','MarkerSize',6)
semilogx(h,OOA(:,2,2),'b--d', 'MarkerFaceColor','b','MarkerSize',6)
semilogx(h,OOA(:,3,2),'g--d', 'MarkerFaceColor','g','MarkerSize',6)
semilogx(h,OOA(:,4,2),'k--d', 'MarkerFaceColor','k','MarkerSize',6)

xlabel('h')
ylabel('Order of Accuracy')
legend('L2-norm rho','L2-norm u', 'L2-norm v', 'L2-norm p',...
    'L1-norm rho', 'L1-norm u', 'L1-norm v', 'L1-norm p')
legend('Location','Best')
grid on


axis([0 16 0 4])











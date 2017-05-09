% Final Project Matlab Plotting 
% Written By: Robert Masti
% This function will read in the data and serve as a template for the different cases
clc; close all; clear all;
num_ghost=3;
NEQ=4;

%% Mesh 1
x9 = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/xc-0.txt');
y9 = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/yc-0.txt');
ni9 = length(x9(1,:)); 
nj9 = length(x9(:,1));

mms9 = zeros(nj9, ni9, NEQ);
num9 = zeros(nj9, ni9, NEQ);
% conserved vars
rhou9 = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/rhou-362.txt');
rhov9 = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/rhov-362.txt');
rhoe9 = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/rhoe-362.txt');
% primvar mms
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/rho_MMS-0.txt');
mms9(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/u_MMS-0.txt');
mms9(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/v_MMS-0.txt');
mms9(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/p_MMS-0.txt');
mms9(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/rho-362.txt');
num9(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/u-362.txt');
num9(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/v-362.txt');
num9(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_1_Flux_1_SS_0_order_2/p-362.txt');
num9(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);


%% Mesh 2
x17 = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/xc-0.txt');
y17 = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/yc-0.txt');
ni17 = length(x17(1,:)); 
nj17 = length(x17(:,1));
mms17 = zeros(nj17, ni17, NEQ);
num17 = zeros(nj17, ni17, NEQ);
% conserved vars
rhou17 = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/rhou-731.txt');
rhov17 = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/rhov-731.txt');
rhoe17 = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/rhoe-731.txt');
% primvar mms
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/rho_MMS-0.txt');
mms17(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/u_MMS-0.txt');
mms17(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/v_MMS-0.txt');
mms17(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/p_MMS-0.txt');
mms17(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/rho-731.txt');
num17(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/u-731.txt');
num17(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/v-731.txt');
num17(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_2_Flux_1_SS_0_order_2/p-731.txt');
num17(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 3
x33 = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/xc-0.txt');
y33 = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/yc-0.txt');
ni33 = length(x33(1,:)); 
nj33 = length(x33(:,1));
mms33 = zeros(nj33, ni33, NEQ);
num33 = zeros(nj33, ni33, NEQ);
% conserved vars
rhou33 = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/rhou-1741.txt');
rhov33 = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/rhov-1741.txt');
rhoe33 = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/rhoe-1741.txt');
% primvar mms
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/rho_MMS-0.txt');
mms33(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/u_MMS-0.txt');
mms33(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/v_MMS-0.txt');
mms33(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/p_MMS-0.txt');
mms33(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/rho-1741.txt');
num33(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/u-1741.txt');
num33(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/v-1741.txt');
num33(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_3_Flux_1_SS_0_order_2/p-1741.txt');
num33(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 4
x65 = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/xc-0.txt');
y65 = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/yc-0.txt');
ni65 = length(x65(1,:)); 
nj65 = length(x65(:,1));
mms65 = zeros(nj65, ni65, NEQ);
num65 = zeros(nj65, ni65, NEQ);
% conserved vars
rhou65 = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/rhou-3627.txt');
rhov65 = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/rhov-3627.txt');
rhoe65 = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/rhoe-3627.txt');
% primvar mms
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/rho_MMS-0.txt');
mms65(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/u_MMS-0.txt');
mms65(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/v_MMS-0.txt');
mms65(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/p_MMS-0.txt');
mms65(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/rho-3627.txt');
num65(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/u-3627.txt');
num65(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/v-3627.txt');
num65(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/p-3627.txt');
num65(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 5
x129 = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/xc-0.txt');
y129 = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/yc-0.txt');
ni129 = length(x129(1,:)); 
nj129 = length(x129(:,1));
mms129 = zeros(nj129, ni129, NEQ);
num129 = zeros(nj129, ni129, NEQ);
% conserved vars
rhou129 = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/rhou-7360.txt');
rhov129 = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/rhov-7360.txt');
rhoe129 = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/rhoe-7360.txt');
% primvar mms
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/rho_MMS-0.txt');
mms129(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/u_MMS-0.txt');
mms129(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/v_MMS-0.txt');
mms129(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/p_MMS-0.txt');
mms129(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/rho-7360.txt');
num129(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/u-7360.txt');
num129(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/v-7360.txt');
num129(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_5_Flux_1_SS_0_order_2/p-7360.txt');
num129(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

%% Mesh 6
x257 = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/xc-0.txt');
y257 = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/yc-0.txt');
ni257 = length(x257(1,:)); 
nj257 = length(x257(:,1));
mms257 = zeros(nj257, ni257, NEQ);
num257 = zeros(nj257, ni257, NEQ);
% conserved vars
rhou257 = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/rhou-16533.txt');
rhov257 = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/rhov-16533.txt');
rhoe257 = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/rhoe-16533.txt');
% primvar mms
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/rho_MMS-0.txt');
mms257(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/u_MMS-0.txt');
mms257(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/v_MMS-0.txt');
mms257(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/p_MMS-0.txt');
mms257(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
% primvar numsoln
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/rho-16533.txt');
num257(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/u-16533.txt');
num257(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/v-16533.txt');
num257(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/p-16533.txt');
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

f3 = figure(3);
f3.PaperUnits = 'inches';
f3.PaperPosition = [0 0 6 5];
h1=subplot(2,2,1);
p1 = get(h1,'pos');
p1(3) = p1(3)+0.075; %right
p1(1) = p1(1)-0.05; %left
p1(4) = p1(4) + 0.05;
p1(2) = p1(2) - 0.05;
set(h1,'pos',p1)
contourf(x257,y257,RelError257(:,:,1) , 50, 'LineStyle', 'none')
title('\rho')
% set(gca,'visible','off')
colorbar
h2=subplot(2,2,2);
p2 = get(h2,'pos');
p2(3) = p2(3)+0.1; %right
p2(1) = p2(1)-0.025; %left
p2(4) = p2(4) + 0.05;
p2(2) = p2(2) - 0.05;
set(h2,'pos',p2)
contourf(x257,y257, RelError257(:,:,2), 50, 'LineStyle', 'none')
title('u')
% set(gca,'visible','off')
colorbar
h3=subplot(2,2,3);
p3 = get(h3,'pos');
p3(3) = p3(3)+0.075; %right
p3(1) = p3(1)-0.05; %left
p3(4) = p3(4) + 0.05;
p3(2) = p3(2) - 0.05;
set(h3,'pos',p3)
contourf(x257,y257,RelError257(:,:,3) , 50, 'LineStyle', 'none')
title('v','FontSize',12)
% set(gca,'visible','off')
colorbar
h4=subplot(2,2,4);
p4 = get(h4,'pos');
p4(3) = p4(3)+0.1; %right
p4(1) = p4(1)-0.025; %left
p4(4) = p4(4) + 0.05;
p4(2) = p4(2) - 0.05;
set(h4,'pos',p4)
contourf(x257,y257, RelError257(:,:,4), 50, 'LineStyle', 'none')
title('P')
% set(gca,'visible','off')
colorbar
 print(f3, 'MMS_mesh6_SB_DE.png', '-dpng', '-r300')

hold off
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

f2 = figure(2);
% f1.PaperUnits = 'inches';
%  f2.PaperPosition = [0 0 8 7];
% L2 norm
semilogx(h,OOA(:,1,1),'r-^', 'MarkerFaceColor','r')
hold on
semilogx(h,OOA(:,2,1),'b-^', 'MarkerFaceColor','b')
semilogx(h,OOA(:,3,1),'g-^', 'MarkerFaceColor','g')
semilogx(h,OOA(:,4,1),'k-^', 'MarkerFaceColor','k')
% Linf norm
semilogx(h,OOA(:,1,2),'r--d', 'MarkerFaceColor','r')
semilogx(h,OOA(:,2,2),'b--d', 'MarkerFaceColor','b')
semilogx(h,OOA(:,3,2),'g--d', 'MarkerFaceColor','g')
semilogx(h,OOA(:,4,2),'k--d', 'MarkerFaceColor','k')

xlabel('h')
ylabel('Order of Accuracy')
leg = legend('L2 \rho','L2 u', 'L2 v', 'L2 P',...
    'L1 \rho', 'L1 u', 'L1 v', 'L1 P');
% set(leg,'position',[0 0 0.2 0.2])
set(leg,'Location','Best')
set(leg,'FontSize',10)

grid on
axis([0 16 0 6])
% print(f2, 'OA_SB.png','-dpng','-r300')
hold off
%% Color Maps

% Plot the numerical Solution

f1 = figure(1);
f1.PaperUnits = 'inches';
f1.PaperPosition = [0 0 6 5];
h1=subplot(2,2,1);
p1 = get(h1,'pos');
p1(3) = p1(3)+0.075; %right
p1(1) = p1(1)-0.05; %left
p1(4) = p1(4) + 0.05;
p1(2) = p1(2) - 0.05;
set(h1,'pos',p1)
contourf(x257,y257,num257(:,:,1) , 50, 'LineStyle', 'none')
title('\rho')
% set(gca,'visible','off')
colorbar
h2=subplot(2,2,2);
p2 = get(h2,'pos');
p2(3) = p2(3)+0.1; %right
p2(1) = p2(1)-0.025; %left
p2(4) = p2(4) + 0.05;
p2(2) = p2(2) - 0.05;
set(h2,'pos',p2)
contourf(x257,y257, num257(:,:,2), 50, 'LineStyle', 'none')
title('u')
% set(gca,'visible','off')
colorbar
h3=subplot(2,2,3);
p3 = get(h3,'pos');
p3(3) = p3(3)+0.075; %right
p3(1) = p3(1)-0.05; %left
p3(4) = p3(4) + 0.05;
p3(2) = p3(2) - 0.05;
set(h3,'pos',p3)
contourf(x257,y257,num257(:,:,3) , 50, 'LineStyle', 'none')
title('v','FontSize',12)
% set(gca,'visible','off')
colorbar
h4=subplot(2,2,4);
p4 = get(h4,'pos');
p4(3) = p4(3)+0.1; %right
p4(1) = p4(1)-0.025; %left
p4(4) = p4(4) + 0.05;
p4(2) = p4(2) - 0.05;
set(h4,'pos',p4)
contourf(x257,y257, num257(:,:,4), 50, 'LineStyle', 'none')
title('P')
% set(gca,'visible','off')
colorbar
% print(f1, 'MMS_mesh6_SB_soln.png', '-dpng', '-r300')





%% Iterative history
% IterEr = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/IterErrHist-3627.txt');
% L2hist = load('Case_1_Mesh_4_Flux_1_SS_0_order_2/L2Hist-3627.txt');
IterEr = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/IterErrHist-16533.txt');
L2hist = load('Case_1_Mesh_6_Flux_1_SS_0_order_2/L2Hist-16533.txt');
for i = 1:4
IterEr(:,i) = IterEr(:,i)/max(max(mms129(:,:,i)));
L2hist(:,i) = L2hist(:,i)/L2hist(1,i);
end


f3 = figure(3);
% f1.PaperUnits = 'inches';
%  f2.PaperPosition = [0 0 8 7];
% L2 norm

semilogy(IterEr(3:end,1),'-r')
hold on
semilogy(IterEr(3:end,2),'-b')
semilogy(IterEr(3:end,3),'-g')
semilogy(IterEr(3:end,4),'-k')
% Linf norm
semilogy(L2hist(3:end,1),'r--')
hold on
semilogy(L2hist(3:end,2),'b--')
semilogy(L2hist(3:end,3),'g--')
semilogy(L2hist(3:end,4),'k--')

xlabel('Iteration')
ylabel('Iterative Residual & Iterative Error')
leg = legend('I.E. \rho','I.E. u', 'I.E. v', 'I.E. P',...
    'L2 \rho', 'L2 u', 'L2 v', 'L2 P');
% set(leg,'position',[0 0 0.2 0.2])
set(leg,'Location','Best')
set(leg,'FontSize',10)

grid on
% print(f3, 'Iter_SB.png','-dpng','-r300')
hold off






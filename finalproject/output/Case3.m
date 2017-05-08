% Final Project Matlab Plotting 
% Written By: Robert Masti
% This function will read in the data and serve as a template for the different cases
clc; close all; clear all;
num_ghost=3;
NEQ=4;

%% Mesh 1 AOA 0

xA0 = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/xc-0.txt');
yA0 = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/yc-0.txt');   
niA0 = length(xA0(1,:)); 
njA0 = length(xA0(:,1));

mesh = 4;
ab = 4*2^mesh+1;
ae = 20*2^mesh;

% Grab the vals
numA0 = zeros(njA0, niA0, NEQ);
% Prim var
temp = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/rho-19967.txt');
numA0(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/u-19967.txt');
numA0(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/v-19967.txt');
numA0(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/p-19967.txt');
numA0(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
%% Mesh 2 AOA 8

xA8 = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/xc-0.txt');
yA8 = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/yc-0.txt');   
niA8 = length(xA8(1,:)); 
njA8 = length(xA8(:,1));

% Grab the vals
numA8 = zeros(njA8, niA8, NEQ);
% Prim var
temp = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/rho-10201.txt');
numA8(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/u-10201.txt');
numA8(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/v-10201.txt');
numA8(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/p-10201.txt');
numA8(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);


%% Plotting
zoom_xmin = -0.15;
zoom_xmax = 0.25;
zoom_ymin = -0.1;
zoom_ymax = 0.1;
clevels = 50;

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
contourf(xA0,yA0,numA0(:,:,1) , clevels,'LineStyle','none'), hold on
plot(xA0(ab:ae,1),yA0(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])
title('\rho')
% set(gca,'visible','off')
colorbar
h2=subplot(2,2,2);
p2 = get(h2,'pos');
p2(3) = p2(3)+0.075; %right
p2(1) = p2(1)-0.0; %left
p2(4) = p2(4) + 0.05;
p2(2) = p2(2) - 0.05;
set(h2,'pos',p2)
contourf(xA0,yA0, numA0(:,:,2), clevels, 'LineStyle','none'), hold on
plot(xA0(ab:ae,1),yA0(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])
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
contourf(xA0,yA0,numA0(:,:,3) , clevels,'LineStyle','none'), hold on
plot(xA0(ab:ae,1),yA0(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])

title('v','FontSize',12)
% set(gca,'visible','off')
colorbar
h4=subplot(2,2,4);
p4 = get(h4,'pos');
p4(3) = p4(3)+0.075; %right
p4(1) = p4(1)-0.0; %left
p4(4) = p4(4) + 0.05;
p4(2) = p4(2) - 0.05;
set(h4,'pos',p4)
contourf(xA0,yA0, numA0(:,:,4), clevels,'LineStyle','none'), hold on
plot(xA0(ab:ae,1),yA0(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])
title('P')
% set(gca,'visible','off')
colorbar
%  print(f1, 'NACA_mesh4_A0_Soln.png', '-dpng', '-r300')

%% Plotting AOA 8
zoom_xmin = -0.15;
zoom_xmax = 0.25;
zoom_ymin = -0.1;
zoom_ymax = 0.2;
clevels = 50;

f2 = figure(2);
f2.PaperUnits = 'inches';
f2.PaperPosition = [0 0 6 5];
h1=subplot(2,2,1);
p1 = get(h1,'pos');
p1(3) = p1(3)+0.075; %right
p1(1) = p1(1)-0.05; %left
p1(4) = p1(4) + 0.05;
p1(2) = p1(2) - 0.05;
set(h1,'pos',p1)
contourf(xA8,yA8,numA8(:,:,1) , clevels,'LineStyle','none'), hold on
plot(xA8(ab:ae,1),yA8(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])
title('\rho')
% set(gca,'visible','off')
colorbar
h2=subplot(2,2,2);
p2 = get(h2,'pos');
p2(3) = p2(3)+0.075; %right
p2(1) = p2(1)-0.0; %left
p2(4) = p2(4) + 0.05;
p2(2) = p2(2) - 0.05;
set(h2,'pos',p2)
contourf(xA8,yA8, numA8(:,:,2), clevels, 'LineStyle','none'), hold on
plot(xA8(ab:ae,1),yA8(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])
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
contourf(xA8,yA8,numA8(:,:,3) , clevels,'LineStyle','none'), hold on
plot(xA8(ab:ae,1),yA8(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])

title('v','FontSize',12)
% set(gca,'visible','off')
colorbar
h4=subplot(2,2,4);
p4 = get(h4,'pos');
p4(3) = p4(3)+0.075; %right
p4(1) = p4(1)-0.0; %left
p4(4) = p4(4) + 0.05;
p4(2) = p4(2) - 0.05;
set(h4,'pos',p4)
contourf(xA8,yA8, numA8(:,:,4), clevels,'LineStyle','none'), hold on
plot(xA8(ab:ae,1),yA8(ab:ae,1),':k')
axis([zoom_xmin zoom_xmax zoom_ymin zoom_ymax])
title('P')
% set(gca,'visible','off')
colorbar
%  print(f2, 'NACA_mesh4_A8_Soln.png', '-dpng', '-r300')


%% Get the Area, lift, and drag forces


% Grab areas to compute components of lift and drag
AjA0 = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/Aj-0.txt');

nj_xA0 = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/n_j_x-0.txt');
nj_yA0 = load('Case_3_Mesh_4_Flux_2_AOA_1_order_1/n_j_y-0.txt');

% Grab A0 forces. 
for j = ab:ae
    F = (numA0(j,1,4)*AjA0(j,1)); % Force Vector
    %mult by -1 to get force on the airfoil not the airfoil on the fluid
    LA0(j,1) = (numA0(j,1,4)*AjA0(j,1))*(-1*nj_yA0(j,1)); 
    DA0(j,1) = (numA0(j,1,4)*AjA0(j,1))*(-1*nj_xA0(j,1));   
end

% Grab areas to compute components of lift and drag
AjA8 = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/Aj-0.txt');

nj_xA8 = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/n_j_x-0.txt');
nj_yA8 = load('Case_3_Mesh_4_Flux_2_AOA_2_order_1/n_j_y-0.txt');

% Grab A8 forces. 
for j = ab:ae
    F = (numA8(j,1,4)*AjA8(j,1)); % Force Vector
    %mult by -1 to get force on the airfoil not the airfoil on the fluid
    LA8(j,1) = (numA8(j,1,4)*AjA8(j,1))*(-1*nj_yA8(j,1)); 
    DA8(j,1) = (numA8(j,1,4)*AjA8(j,1))*(-1*nj_xA8(j,1));   
end

output = ['AOA 0: D = ', num2str(sum(DA0)), ' , L = ', num2str(sum(LA0))];
output2 = ['AOA 8: D = ', num2str(sum(DA8)), ' , L = ', num2str(sum(LA8))];

disp(output), disp(output2)

% Final Project Matlab Plotting 
% Written By: Robert Masti
% This function will read in the data and serve as a template for the different cases
clc; close all; clear all;
num_ghost=3;
NEQ=4;

%% Mesh 1

x = load('Case_2_Mesh_4_Flux_1_order_1/xc-0.txt');
y = load('Case_2_Mesh_4_Flux_1_order_1/yc-0.txt');
ni = length(x(1,:)); 
nj = length(x(:,1));

% Grab the vals
mms = zeros(nj, ni, NEQ);
num = zeros(nj, ni, NEQ);
% Prim var
temp = load('Case_2_Mesh_4_Flux_1_order_1/rho-1437.txt');
num(:,:,1) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_2_Mesh_4_Flux_1_order_1/u-1437.txt');
num(:,:,2) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_2_Mesh_4_Flux_1_order_1/v-1437.txt');
num(:,:,3) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);
temp = load('Case_2_Mesh_4_Flux_1_order_1/p-1437.txt');
num(:,:,4) = temp(num_ghost+1:end-num_ghost,num_ghost+1:end-num_ghost);

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
contourf(x,y,num(:,:,1) , 50, 'LineStyle', 'none')
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
contourf(x,y, num(:,:,2), 50, 'LineStyle', 'none')
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
contourf(x,y,num(:,:,3) , 50, 'LineStyle', 'none')
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
contourf(x,y, num(:,:,4), 50, 'LineStyle', 'none')
title('P')
% set(gca,'visible','off')
colorbar
% print(f1, 'Inlet_mesh4_Soln.png', '-dpng', '-r300')


%% Get the Average Pressure at the outlet

A = load('Case_2_Mesh_4_Flux_1_order_1/Ai-0.txt');

for n = 1:length(x(nj,:))
    F(n,1) = num(nj,n,4)*A(end,n);
end
avgF = sum(F)/nj
avgP = avgF/(sum(A(end,:)))



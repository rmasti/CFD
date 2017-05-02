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
rhouc2d9 = load('Case_1_Mesh_1_Flux_1_order_2/rhou-79.txt');
rhovc2d9 = load('Case_1_Mesh_1_Flux_1_order_2/rhov-79.txt');
rhoec2d9 = load('Case_1_Mesh_1_Flux_1_order_2/rhoe-79.txt');
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





% contourf(xc,yc,rhoc)
% colorbar

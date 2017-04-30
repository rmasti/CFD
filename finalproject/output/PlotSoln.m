% Final Project Matlab Plotting 
% Written By: Robert Masti
% This function will read in the data and serve as a template for the different cases
clc; close all; clear all;
num_ghost=3;
xc = load('Case_1_Mesh_3_Flux_1_order_2/xc-0.txt');
yc = load('Case_1_Mesh_3_Flux_1_order_2/yc-0.txt');
rhocg = load('Case_1_Mesh_3_Flux_1_order_2/rho-100.txt');

rhoc = rhocg(num_ghost+1:end-num_ghost, num_ghost+1:end-num_ghost);

contourf(xc,yc,rhoc)
colorbar

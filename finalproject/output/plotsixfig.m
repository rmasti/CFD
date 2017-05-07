% Plotting script for HLLC euler equation
% Written By: Robert Masti
clc, clear; close all;

% Set the number of ghost cells
ng = 3;

width = 3;     % Width in inches 
height = 3;    % Height in inches 
alw = 0.75;    % AxesLineWidth 
fsz = 11;      % Fontsize 
lw = 1.5;      % LineWidth 
msz = 8;       % MarkerSize 
% load the data
load('hllc_euler.mat');

% Plot the time evolution of density using a subplot
f1 = figure(1);
f1.PaperUnits = 'inches'
f1.PaperPosition = [0 0 6 6]
h1=subplot(3,2,1);
p1 = get(h1,'pos')
p1(3) = p1(3)+0.095; %right
p1(1) = p1(1)-0.075;
p1(4) = p1(4) + 0.075;
p1(2) = p1(2) - 0.05;
set(h1,'pos',p1)
contourf(history(40).U(1+ng:end-ng,1+ng:end-ng,3), 50, 'LineStyle', 'none')
set(gca,'visible','off')
colormap('hot')

% colorbar('Position', [p6(3)+p5(3)+0.15  p6(2)  0.025  p6(4)+p5(4)*2.4])
% colorbar('Position', [p1(3)+p1(3)+0.15  -0.25*p1(2)  0.025  p1(4)+p1(4)*2.4])
h2=subplot(3,2,3);
p2 = get(h2,'pos')
p2(3) = p2(3)+0.095; %right
p2(1) = p2(1)-0.075; %left
p2(4) = p2(4) + 0.075;
p2(2) = p2(2) - 0.05;


set(h2,'pos',p2)
contourf(history(60).U(1+ng:end-ng,1+ng:end-ng,3), 50, 'LineStyle', 'none')
set(gca,'visible','off')
colormap('hot')



h3=subplot(3,2,5);
p3 = get(h3,'pos')
p3(3) = p3(3)+0.095;
p3(1) = p3(1)-0.075;
p3(4) = p3(4) + 0.075;
p3(2) = p3(2) - 0.05;
set(h3,'pos',p3)
contourf(history(80).U(1+ng:end-ng,1+ng:end-ng,3), 50, 'LineStyle', 'none')
set(gca,'visible','off')
colormap('hot')

h4=subplot(3,2,2);
p4 = get(h4,'pos')
p4(3) = p4(3)+0.075; %right
p4(1) = p4(1)-0.075;
p4(4) = p4(4) + 0.075;
p4(2) = p4(2) - 0.05;

set(h4,'pos',p4)
contourf(history(120).U(1+ng:end-ng,1+ng:end-ng,3), 50, 'LineStyle', 'none')
set(gca,'visible','off')
colormap('hot')


h5=subplot(3,2,4);
p5 = get(h5,'pos')
p5(3) = p5(3)+0.075; %right
p5(1) = p5(1)-0.075; %left
p5(4) = p5(4) + 0.075;
p5(2) = p5(2) - 0.05;

set(h5,'pos',p5)
contourf(history(140).U(1+ng:end-ng,1+ng:end-ng,3), 50, 'LineStyle', 'none')
set(gca,'visible','off')
colormap('hot')



h6=subplot(3,2,6);
p6 = get(h6,'pos')
p6(3) = p6(3)+0.075;
p6(1) = p6(1)-0.075;
p6(4) = p6(4) + 0.075;
p6(2) = p6(2) - 0.05;

set(h6,'pos',p6)
contourf(history(160).U(1+ng:end-ng,1+ng:end-ng,3), 50, 'LineStyle', 'none')
set(gca,'visible','off')
colormap('hot')

colorbar('Position', [p6(3)+p5(3)+0.095  p6(2)  0.025  p6(4)+p5(4)*2.0])
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
title =text(0.5, 0.995,'\bf y-momentum ','HorizontalAlignment','center','VerticalAlignment', 'top');
x = text(0.5, 0.05,'x','HorizontalAlignment','center','VerticalAlignment', 'top');
y = text(0.0, 0.5,'y','HorizontalAlignment','center','VerticalAlignment', 'top');
set(y, 'Rotation', 90)
set(gca,'visible','off')
% x.FontSize(12), title.FontSize(14), y.FontSize(12)

print(f1, 'hllc_time_evolution2_ymtm.png', '-dpng', '-r300')






















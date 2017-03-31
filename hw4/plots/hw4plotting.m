% Plotting Script Hw2 CFD
% Written By: Robert Masti
% 02/27/2017
clc, clear, close all
% IMPORT SOLUTION HAS THE FOLLOWING VARIABLES
% x A rho u v p M U1 U2 U3 U4
% 1 2  3  4 5 6 7 8  9  10 11
%% Effect of Kappa2

% Import solutions
fourthsoln = importsolution('kappa2/fourth/sb/q1Dnozzle.dat');
thirdsoln = importsolution('kappa2/third/sb/q1Dnozzle.dat');
halfsoln = importsolution('kappa2/half/sb/q1Dnozzle.dat');

f1 = figure();
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.75, 9.125], ...
'PaperUnits', 'Inches', 'PaperSize', [7.75, 9.125])
subplot(3,1,1)
plot(thirdsoln(:,1),thirdsoln(:,8),'r',...
    halfsoln(:,1),halfsoln(:,8),'--b')
axis([0.4 0.6 0.16 0.26])
legend('\kappa_2 = 1/3', '\kappa_2 = 1/2')
ylabel('density'), grid on
title('Effect of Kappa2 on Shock Capturing (256 cells)')
subplot(3,1,2)
plot(thirdsoln(:,1),thirdsoln(:,9),'r',...
    halfsoln(:,1),halfsoln(:,9),'--b')
legend('\kappa_2 = 1/3', '\kappa_2 = 1/2')
ylabel('mtm'), grid on
axis([0.4 0.6 140 210])
subplot(3,1,3)
plot(thirdsoln(:,1),thirdsoln(:,11),'r',...
    halfsoln(:,1),halfsoln(:,11),'--b')
legend('\kappa_2 = 1/3', '\kappa_2 = 1/2')
ylabel('energy'), grid on
xlabel('x')
axis([0.44 0.56 0.8e5 1.6e5])
filename = 'plots/kappa2.png';
saveas(f1,filename)

%% Effect of Kappa4

% Import solutions
sixtyfourthsoln = importsolution('kappa4/64th/sb/q1Dnozzle.dat');
fortyeighthsoln = importsolution('kappa4/48th/sb/q1Dnozzle.dat');
thirtysecondsoln = importsolution('kappa4/32nd/sb/q1Dnozzle.dat');

f2 = figure();
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.75, 9.125], ...
'PaperUnits', 'Inches', 'PaperSize', [7.75, 9.125])
subplot(3,1,1)
plot(sixtyfourthsoln(:,1),sixtyfourthsoln(:,8),'r',...
    fortyeighthsoln(:,1),fortyeighthsoln(:,8),'--b', ...
    thirtysecondsoln(:,1), thirtysecondsoln(:,8), '-.k')
axis([0.54 0.62 0.66 0.68])
legend('\kappa_4 = 1/64', '\kappa_4 = 1/48', '\kappa_4 = 1/32')
ylabel('density'), grid on
title('Effect of Kappa4 on Shock Capturing (256 cells)')
subplot(3,1,2)
plot(sixtyfourthsoln(:,1),sixtyfourthsoln(:,9),'r',...
    fortyeighthsoln(:,1),fortyeighthsoln(:,9),'--b', ...
    thirtysecondsoln(:,1), thirtysecondsoln(:,9), '-.k')
legend('\kappa_4 = 1/64', '\kappa_4 = 1/48', '\kappa_4 = 1/32')
ylabel('mtm'), grid on
axis([0.54 0.6 138 150])
subplot(3,1,3)
plot(sixtyfourthsoln(:,1),sixtyfourthsoln(:,11),'r',...
    fortyeighthsoln(:,1),fortyeighthsoln(:,11),'--b', ...
    thirtysecondsoln(:,1), thirtysecondsoln(:,11), '-.k')
legend('\kappa_4 = 1/64', '\kappa_4 = 1/48', '\kappa_4 = 1/32')
ylabel('energy'), grid on
xlabel('x')
axis([0.52 0.62 2.85e5 3.0e5])
filename = 'plots/kappa4.png';
saveas(f2,filename)

%% Order of Accuracy

% Import soln and exact soln

g32 = importsolution('meshrefine/32/ss/q1Dnozzle.dat');
g32e = importsolution('meshrefine/32/ss/exactsol.dat');
DE32(1) = sqrt(sum((g32(:,8)-g32e(:,8)).^2)/length(g32(:,8)));
DE32(2) = sqrt(sum((g32(:,9)-g32e(:,9)).^2)/length(g32(:,9)));
DE32(3) = sqrt(sum((g32(:,11)-g32e(:,11)).^2)/length(g32(:,11)));

g64 = importsolution('meshrefine/64/ss/q1Dnozzle.dat');
g64e = importsolution('meshrefine/64/ss/exactsol.dat');
DE64(1) = sqrt(sum((g64(:,8)-g64e(:,8)).^2)/length(g64(:,8)));
DE64(2) = sqrt(sum((g64(:,9)-g64e(:,9)).^2)/length(g64(:,9)));
DE64(3) = sqrt(sum((g64(:,11)-g64e(:,11)).^2)/length(g64(:,11)));

g128 = importsolution('meshrefine/128/ss/q1Dnozzle.dat');
g128e = importsolution('meshrefine/128/ss/exactsol.dat');
DE128(1) = sqrt(sum((g128(:,8)-g128e(:,8)).^2)/length(g128(:,8)));
DE128(2) = sqrt(sum((g128(:,9)-g128e(:,9)).^2)/length(g128(:,9)));
DE128(3) = sqrt(sum((g128(:,11)-g128e(:,11)).^2)/length(g128(:,11)));

g256 = importsolution('meshrefine/256/ss/q1Dnozzle.dat');
g256e = importsolution('meshrefine/256/ss/exactsol.dat');
DE256(1) = sqrt(sum((g256(:,8)-g256e(:,8)).^2)/length(g256(:,8)));
DE256(2) = sqrt(sum((g256(:,9)-g256e(:,9)).^2)/length(g256(:,9)));
DE256(3) = sqrt(sum((g256(:,11)-g256e(:,11)).^2)/length(g256(:,11)));

g512 = importsolution('meshrefine/512/ss/q1Dnozzle.dat');
g512e = importsolution('meshrefine/512/ss/exactsol.dat');
DE512(1) = sqrt(sum((g512(:,8)-g512e(:,8)).^2)/length(g512(:,8)));
DE512(2) = sqrt(sum((g512(:,9)-g512e(:,9)).^2)/length(g512(:,9)));
DE512(3) = sqrt(sum((g512(:,11)-g512e(:,11)).^2)/length(g512(:,11)));


for i=1:3
    
    p(i,1) = log(DE256(i)/DE512(i))/log(2);
    p(i,2) = log(DE128(i)/DE256(i))/log(2);
    p(i,3) = log(DE64(i)/DE128(i))/log(2);
    p(i,4) = log(DE32(i)/DE64(i))/log(2);
    
end
h = [ 1 4 8 16 ];

f3 = figure();
semilogx(h, p(1,:) , '-*r',...
    h, p(2,:) , '-db',...
    h, p(2,:) , '-sk')
legend('L2 density', 'L2 mtm', 'L2 energy'), xlabel('h')
ylabel('p'), title('Observed Order of Accuracy')
grid on, axis([0 20 0 4])
filename = 'plots/order.png';
saveas(f3,filename)
%% Demonstrate Iterative Convergence
% Iteration Time L2norm1 L2norm2 L2norm3 L2norm4
% import data
h32 = importsolution('meshrefine/32/ss/q1Dhistory.dat');
h64 = importsolution('meshrefine/64/ss/q1Dhistory.dat');
h128 = importsolution('meshrefine/128/ss/q1Dhistory.dat');
h256 = importsolution('meshrefine/256/ss/q1Dhistory.dat');
h256sb = importsolution('meshrefine/256/sb/q1Dhistory.dat');
h512 = importsolution('meshrefine/512/ss/q1Dhistory.dat');
f4 = figure();
subplot(2,1,1)
semilogy(h256(:,1), h256(:,3), '-r', ...
    h256(:,1), h256(:,4), '--b', ...
    h256(:,1), h256(:,6), '-.k')
legend('L2 density', 'L2 mtm', 'L2 energy')
xlabel('iteration'), grid on
ylabel('residual')
title('Iterative Convergence Supersonic Outflow (256 cells)')
subplot(2,1,2)
semilogy(h256sb(:,1), h256sb(:,3), '-r', ...
    h256sb(:,1), h256sb(:,4), '--b', ...
    h256sb(:,1), h256sb(:,6), '-.k')
legend('L2 density', 'L2 mtm', 'L2 energy')
xlabel('iteration'), grid on
ylabel('residual')
title('Iterative Convergence Subsonic Outflow (256 cells)')
filename = 'plots/conv.png';
saveas(f4,filename)

%% Grid Effect on soltn

% g32 g64 g128 g256 g512

f5 = figure();
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.75, 9.125], ...
'PaperUnits', 'Inches', 'PaperSize', [7.75, 9.125])
subplot(3,1,1)
plot(g512e(:,1), g512e(:,8), 'r'), hold on
plot(g32(:,1), g32(:,8), 'db',...
    g128(:,1), g128(:,8), '.k')
legend('Exact','32', '128'), grid on
ylabel('density'), title('Grid Size Effect on Supersonic Outflow Solution')
hold off
subplot(3,1,2)
plot(g512e(:,1), g512e(:,9), 'r'), hold on
plot(g32(:,1), g32(:,9), 'db',...
    g128(:,1), g128(:,9), '.k')
legend('Exact','32', '128')
ylabel('mtm'), grid on
hold off
subplot(3,1,3)
plot(g512e(:,1), g512e(:,11), 'r'), hold on
plot(g32(:,1), g32(:,11), 'db',...
    g128(:,1), g128(:,11), '.k')
legend('Exact','32', '128')
ylabel('mtm'), grid on, xlabel('x')
hold off
filename = 'plots/gridss.png';
saveas(f5,filename)

% import subsonic solution
g32sb = importsolution('meshrefine/32/sb/q1Dnozzle.dat');
g128sb = importsolution('meshrefine/128/sb/q1Dnozzle.dat');
g512sb = importsolution('meshrefine/512/sb/q1Dnozzle.dat');

f6 = figure();
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.75, 9.125], ...
'PaperUnits', 'Inches', 'PaperSize', [7.75, 9.125])
subplot(3,1,1)
plot(g512sb(:,1), g512sb(:,8), 'r'), hold on
plot(g32sb(:,1), g32sb(:,8), 'db',...
    g128sb(:,1), g128sb(:,8), '.k')
legend('512','32', '128'), grid on
ylabel('density'), title('Grid Size Effect on Subsonic Outflow Solution')
hold off
subplot(3,1,2)
plot(g512sb(:,1), g512sb(:,9), 'r'), hold on
plot(g32sb(:,1), g32sb(:,9), 'db',...
    g128sb(:,1), g128sb(:,9), '.k')
legend('512','32', '128')
ylabel('mtm'), grid on
hold off
subplot(3,1,3)
plot(g512sb(:,1), g512sb(:,11), 'r'), hold on
plot(g32sb(:,1), g32sb(:,11), 'db',...
    g128sb(:,1), g128sb(:,11), '.k')
legend('512','32', '128')
ylabel('mtm'), grid on, xlabel('x')
hold off
filename = 'plots/gridsb.png';
saveas(f6,filename)




















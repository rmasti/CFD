# CFD
computational fluid dynamics
This repo will contain all c++ files involved in computational fluid dynamics course offered by virginia tech through the AOE department. 

Professor: Dr. Roy

Creator: R. Masti

Date 02/13/2017

Structure details will be added here as the code is developed. 
NOTE when compiling one must use c++ 11 example:

[user directory]$ g++ -std=c++11x "all cpp files" -o "executable"

The folder hw2 contains the code used to model an quasi-linear axisymmetric 1D flow through a nozzle it incorporates jameson damping and using a simple Finite Volume Method scheme with cell averaged values. The code also uses a Newton-Raphson method to find the exact solution by computing the mach number through the isentropic relation. The code also has two version v1 and v2 and v2 is the newer version that should be adapted becuase of its user friendly approach. They both use a unit testing framework call cppCatch, and only v2 uses a module called Eigen which is used for having matrices and doing manipulation similar to what is possible in MATLAB. Also there is a MakeFile located in the main folder and can be used for easy compiling purposes. To compile tests just run make tests, and run make for the main executable.

NOTE to save space the shared package called Eigen was put in the top level directory which is both used by hw4 an hw2 which is required. There is an issue also with the output handling and making sure that the directory in which it wants to output in actually exists. So to counter this the directors data/ss and data/sb have a README file which allows them to not be empty and therefore will be committed to the repo

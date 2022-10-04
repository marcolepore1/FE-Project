function Parameters=FFTparameters(M,dz_x1, flag)
%
% Function that computes the parameters of the Fast Fourier transform grid
% given in input one between dz or x1
%
% INPUT
% 2^M:                  number of discretizations
% dz:                   interval of 
% flag:                 flag==1 if we pass dz otherwise we pass x1
%
% OUTPUT
% Parameters:           struct with FFT parameters
%
%

N=2^M;
if flag == 1
    dz = dz_x1;
    dx = 2*pi/(N * dz);
    x1 = - 0.5*(N-1)*dx;
else 
    x1 = dz_x1;
    dx = - 2*x1 / (N-1);
    dz = 2*pi/(N * dx);
end
    
x  = [x1:dx:-x1];
z1 = - 0.5*(N-1)*dz;
z  = [z1:dz:-z1];
Parameters.x1 = z1;
Parameters.xn = -z1;
Parameters.dx = dz;
Parameters.z1 = x1;
Parameters.zn = -x1;
Parameters.dz = dx;
Parameters.M  = M;
    
end
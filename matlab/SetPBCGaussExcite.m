function FDTD = SetPBCGaussExcite(FDTD,f0,fc)
% function FDTD = SetPBCGaussExcite(FDTD,f0,fc);
%
% f0 : center frequency
% fc : 20dB cutoff frequency --> bandwidth is 2*fc
%
% see also SetSinusExcite SetCustomExcite
%
% e.g FDTD = SetPBCGaussExcite(FDTD,1e9,1e8);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.Excitation.ATTRIBUTE.Type=4;
FDTD.Excitation.ATTRIBUTE.f0=f0;
FDTD.Excitation.ATTRIBUTE.fc=fc;
FDTD.ATTRIBUTE.f_max=f0+fc;

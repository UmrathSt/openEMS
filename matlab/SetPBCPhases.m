function FDTD = SetPBCPhases(FDTD, phases)
% FDTD = SetPBCPhases(FDTD, phases)
%
% phases = [phase_x, phase_y, phase_z]
%
% phase_x is the phase which occurs upon propagation
% of the incident wave through the studied structure
% phase_x = L_x * k_x
% with the wave-vector component k_x 
% openEMS matlab interface
% -----------------------
% author: Stefan Umrath

if (numel(phases)~=3)
    error('openEMS:SetPBCPhases','wrong number of phasefactors');
end

if isnumeric(phases)
    FDTD.PBCPhases.ATTRIBUTE.pbc_phase_x=phases(1);
    FDTD.PBCPhases.ATTRIBUTE.pbc_phase_y=phases(2);
    FDTD.PBCPhases.ATTRIBUTE.pbc_phase_z=phases(3);
elseif iscell(phases)
    FDTD.PBCPhases.ATTRIBUTE.pbc_phase_x=phases{1};
    FDTD.PBCPhases.ATTRIBUTE.pbc_phase_y=phases{2};
    FDTD.PBCPhases.ATTRIBUTE.pbc_phase_z=phases{3};
else
    error('openEMS:SetPBCPhases','unknown type of phasefactors');
end

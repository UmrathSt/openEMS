function [CSX,port] = AddPBCWaveGuidePort( CSX, prio, portnr, start, stop, dir, E_WG_funcsin,E_WG_funccos, H_WG_funcsin,H_WG_funccos, exc_amp, exc_type, varargin )
% function [CSX,port] = AddPBCWaveGuidePort( CSX, prio, portnr, start, stop, dir, E_WG_func, H_WG_func, exc_amp, exc_type, varargin )
% 
% Create a PBC waveguide port, including an optional excitation and probes
% 
% Note: - The excitation will be located at the start position in the given direction
%       - The voltage and current probes at the stop position in the given direction
%
% parameter:
%   CSX:        complete CSX structure (must contain a mesh)
%   prio:       priority of primitives
%   start:      start coordinates of waveguide port box
%   stop:       stop  coordinates of waveguide port box
%   dir:        direction of port (0/1/2 or 'x'/'y'/'z'-direction)
%   E_WG_funcsin:  electric field mode profile function as a string for sin(t)
%   E_WG_funccos:  electric field mode profile function as a string for cos(t)
%   H_WG_funcsin:  magnetic field mode profile function as a string for sin(t)
%   H_WG_funccos:  magnetic field mode profile function as a string for cos(t)
%   exc_amp:    excitation amplitude (set 0 to be passive)
%   exc_type    type of the excitation 0: -> E-soft-excite, 2-> H-soft-excite
%
% optional (key/values):
%   varargin:   optional additional excitations options, see also AddExcitation
%   'PortNamePrefix': a prefix to the port name
%
% output:
%   CSX:        modified CSX structure
%   port:       port structure to use with calcPort
%
% example:
%   % create a TE11 circular waveguide mode, using cylindircal coordinates
%   p11 = 1.841;
%
%   start=[mesh.r(1)   mesh.a(1)   0  ];
%   stop =[mesh.r(end) mesh.a(end) 100];
%   [CSX, port{1}] = AddPBCWaveGuidePort(CSX, 0, 1, start, stop, 2, func_E, func_H, 1);
%
% openEMS matlab interface
% -----------------------
% (c) 2013 Thorsten Liebig (thorsten.liebig@gmx.de)
%
% See also InitCSX, AddExcitation, calcWGPort, calcPort

%check mesh
if ~isfield(CSX, 'RectilinearGrid')
    error('mesh needs to be defined! Use DefineRectGrid() first!');
end

dir = DirChar2Int(dir);

port.type='PBCWaveGuide';
port.nr=portnr;
port.dir = dir;
port.drawingunit = CSX.RectilinearGrid.ATTRIBUTE.DeltaUnit;

PortNamePrefix = 'PBC';

varargin_tmp  = varargin;
for n=1:2:numel(varargin_tmp)
    if strcmpi('PortNamePrefix',varargin_tmp{n})
        PortNamePrefix = varargin_tmp{n+1};
        varargin([n n+1]) = [];
    end
end

% matlab adressing
dir = dir + 1;
dir_sign = sign(stop(dir) - start(dir));
if (dir_sign==0)
    dir_sign = 1;
end

port.direction = dir_sign;


port.excite = 0;
if (exc_amp~=0)
    if (start(dir)==stop(dir))
        error 'if waveguide port is to be excited, the length in propagation direction must not be zero'
    end
    e_start = start;
    e_stop = stop;
    port.excite = 1;
    port.excitepos = e_start(dir);
    e_vec = [1 1 1]*exc_amp;
    e_vec(dir) = 0;
    exc_nameE = [PortNamePrefix '_port_exciteE_' num2str(portnr)];
    exc_nameH = [PortNamePrefix '_port_exciteH_' num2str(portnr)];
    if ~(exc_type==0 || exc_type==2)
        error('The excitation type for PBC excitation must be either 0 (softE) or 2 (softH) but got: %i', exc_type);
    end
    fprintf('I am now adding current an voltage excitations for the PBC');
    CSX = AddPBCExcitation( CSX, exc_nameE, 0, e_vec, varargin{:});
    CSX = AddPBCExcitation( CSX, exc_nameH, 2, e_vec, varargin{:});
    CSX = SetPBCExcitationWeight(CSX, exc_nameE, E_WG_funcsin, E_WG_funccos);
    CSX = SetPBCExcitationWeight(CSX, exc_nameH, H_WG_funcsin, H_WG_funccos);
    
	CSX = AddBox( CSX, exc_nameE, prio, e_start, e_stop);
    CSX = AddBox( CSX, exc_nameH, prio, e_start, e_stop);
end

% voltage/current planes
m_start = start;
m_stop = stop;
m_start(dir) = stop(dir);

port.measplanepos = m_start(dir);
port.U_filename = [PortNamePrefix 'port_ut' int2str(portnr)];
CSX = AddProbe(CSX, port.U_filename, 10, 'ModeFunction', E_WG_funcsin);
CSX = AddBox(CSX, port.U_filename, 0 ,m_start, m_stop);

port.I_filename = [PortNamePrefix 'port_it' int2str(portnr)];
CSX = AddProbe(CSX, port.I_filename, 11, 'ModeFunction', H_WG_funcsin, 'weight', dir_sign);
CSX = AddBox(CSX, port.I_filename, 0 ,m_start, m_stop);
end

/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*   This module has been implemented by Stefan Umrath at the
*   German Aerospace Center (DLR e. V., mail: Stefan.Umrath@dlr.de )
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENGINE_EXT_PBC_H
#define ENGINE_EXT_PBC_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

class Operator_Ext_Pbc;

class Engine_Ext_Pbc : public Engine_Extension
{
    friend class Engine;
    friend class Operator_Ext_Pbc;
public:
    Engine_Ext_Pbc(Operator_Ext_Pbc* op_ext);
    virtual ~Engine_Ext_Pbc();

    virtual void SetNumberOfThreads(int nrThread);

    virtual void DoPreVoltageUpdates() {Engine_Ext_Pbc::DoPreVoltageUpdates(0);}
    virtual void DoPreVoltageUpdates(int threadID);

    virtual void DoPreCurrentUpdates() {Engine_Ext_Pbc::DoPreCurrentUpdates(0);}
    virtual void DoPreCurrentUpdates(int threadID);


    virtual void DoPostVoltageUpdates() {Engine_Ext_Pbc::DoPostVoltageUpdates(0);}
    virtual void DoPostVoltageUpdates(int threadID);

    virtual void DoPostCurrentUpdates() {Engine_Ext_Pbc::DoPostCurrentUpdates(0);};
    virtual void DoPostCurrentUpdates(int threadID);

    void DoVoltPhaseUpdates();
    void DoCurrPhaseUpdates();
    void Apply2Voltages();

    void Apply2Current();






protected:
    Operator_Ext_Pbc* m_Op_Pbc;
    inline bool IsActive() {if (m_Eng->GetNumberOfTimesteps()<m_start_TS) return false; return true;};
    unsigned int m_start_TS;
    int m_ny;
    int m_nyP,m_nyPP;
    unsigned int m_LineNr;
    int m_LineNr_Shift;
    unsigned int m_numLines[3];
    bool direction_is_pbc[6];
    FDTD_FLOAT *pbc_phase;
    bool shift[3];
    unsigned int pos[3];

    vector<unsigned int> m_start;
    vector<unsigned int> m_numX;

    FDTD_FLOAT**** volt_im; // imaginary part of the complex voltage
    FDTD_FLOAT**** curr_im; // imaginary part of the complex current
    FDTD_FLOAT sin_kx;
    FDTD_FLOAT sin_ky;
    FDTD_FLOAT cos_kx;
    FDTD_FLOAT cos_ky;
    FDTD_FLOAT curr_im_outside[3] = {0};
    FDTD_FLOAT curr_outside[3] = {0};
    FDTD_FLOAT volt_im_outside[3] = {0};
    FDTD_FLOAT volt_outside[3] = {0};
    unsigned int maxX, maxY;
};

#endif // ENGINE_EXT_PBC_H

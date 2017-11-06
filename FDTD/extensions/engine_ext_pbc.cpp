/* engine extension for periodic boundary conditions
*/
#include "engine_ext_pbc.h"
#include "operator_ext_pbc.h"
#include "operator_extension.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "operator_ext_excitation.h"

#include "FDTD/engine.h"

Engine_Ext_Pbc::Engine_Ext_Pbc(Operator_Ext_Pbc* op_ext) : Engine_Extension(op_ext)
{
    m_Op_Pbc = op_ext;
    m_numLines[0] = m_Op_Pbc->m_numLines[0];
    m_numLines[1] = m_Op_Pbc->m_numLines[1];
    m_numLines[2] = m_Op_Pbc->m_numLines[2];
    volt_im = Create_N_3DArray<FDTD_FLOAT>(m_numLines); // imaginary parts of
    curr_im = Create_N_3DArray<FDTD_FLOAT>(m_numLines); // voltage/current as (3, Nx, Ny, Nz) array of floats
    SetNumberOfThreads(1);
    for (unsigned int i = 0; i<3; ++i){
        direction_is_pbc[i] = m_Op_Pbc->m_Op->dir_is_pbc[i];
    }

}


Engine_Ext_Pbc::~Engine_Ext_Pbc()
{
    Delete_N_3DArray(volt_im,m_numLines);
    volt_im=NULL;
    Delete_N_3DArray(curr_im,m_numLines);
    curr_im=NULL;
}

void Engine_Ext_Pbc::SetNumberOfThreads(int nrThread)
{
    Engine_Extension::SetNumberOfThreads(nrThread);

    m_numX = AssignJobs2Threads(m_numLines[0],m_NrThreads,false);
    m_start.resize(m_NrThreads,0);
    m_start.at(0)=0;
    for (size_t n=1; n<m_numX.size(); ++n)
        m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}

void Engine_Ext_Pbc::DoPostVoltageUpdates(int threadID){
    unsigned int pos[3];
    bool shift[3];

    for (pos[0]=0; pos[0]<m_numLines[0]-1; ++pos[0])
    {
        for (pos[1]=0; pos[1]<m_numLines[1]-1; ++pos[1])
        {
            for (pos[2]=0; pos[2]<m_numLines[2]-1; ++pos[2])
            {
                //do the updates here
                //for x
                volt_im[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(0,pos[0],pos[1],pos[2]);
                volt_im[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(0,pos[0],pos[1],pos[2]) * (curr_im[2][pos[0]][pos[1]][pos[2]] - curr_im[2][pos[0]][pos[1]+1][pos[2]] - curr_im[1][pos[0]][pos[1]][pos[2]] + curr_im[1][pos[0]][pos[1]][pos[2]+1]);

                //for y
                volt_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(1,pos[0],pos[1],pos[2]);
                volt_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(1,pos[0],pos[1],pos[2]) * (curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im[0][pos[0]][pos[1]][pos[2]+1] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im[2][pos[0]+1][pos[1]][pos[2]]);

                //for z
                volt_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(2,pos[0],pos[1],pos[2]);
                volt_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(2,pos[0],pos[1],pos[2]) * (curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im[1][pos[0]+1][pos[1]][pos[2]]- curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im[0][pos[0]][pos[1]+1][pos[2]]);

            }
        }
    }
};

void Engine_Ext_Pbc::DoPostCurrentUpdates(int threadID){
    unsigned int pos[3];
    bool shift[3];

    for (pos[0]=0; pos[0]<m_numLines[0]-1; ++pos[0])
    {
        for (pos[1]=0; pos[1]<m_numLines[1]-1; ++pos[1])
        {
            for (pos[2]=0; pos[2]<m_numLines[2]-1; ++pos[2])
            {
                //do the updates here
                //for x
                curr_im[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(0,pos[0],pos[1],pos[2]);
                curr_im[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(0,pos[0],pos[1],pos[2]) * (volt_im[2][pos[0]][pos[1]][pos[2]] - volt_im[2][pos[0]][pos[1]+1][pos[2]] - volt_im[1][pos[0]][pos[1]][pos[2]] + volt_im[1][pos[0]][pos[1]][pos[2]+1]);

                //for y
                curr_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(1,pos[0],pos[1],pos[2]);
                curr_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(1,pos[0],pos[1],pos[2]) * (volt_im[0][pos[0]][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]+1] - volt_im[2][pos[0]][pos[1]][pos[2]] + volt_im[2][pos[0]+1][pos[1]][pos[2]]);

                //for z
                curr_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(2,pos[0],pos[1],pos[2]);
                curr_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(2,pos[0],pos[1],pos[2]) * (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im[1][pos[0]+1][pos[1]][pos[2]]- volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im[0][pos[0]][pos[1]+1][pos[2]]);

            }
        }
    }
};
void Engine_Ext_Pbc::Apply_Phases_to_dir(unsigned int dir){
    m_ny   = dir;
    m_nyP  = (dir+1)%3;
    m_nyPP = (dir+2)%3;
    unsigned int posL[3];
    unsigned int posR[3];
    unsigned int dir_lines[2] = {0, m_numLines[dir]-1};
    FDTD_FLOAT tmp_Iim[3];
    FDTD_FLOAT tmp_Ire[3];
    FDTD_FLOAT tmp_Uim[3];
    FDTD_FLOAT tmp_Ure[3];
    FDTD_FLOAT sinus = sin(m_Op_Pbc->m_Op->m_k_PBC[dir]);
    FDTD_FLOAT cosinus = cos(m_Op_Pbc->m_Op->m_k_PBC[dir]);

    posL[m_ny] = dir_lines[0];
    posR[m_ny] = dir_lines[1];
    for (posL[m_nyP]=0; posL[m_nyP]<m_numLines[m_nyP]-1; ++posL[m_nyP]){
        for (posL[m_nyPP]=0; posL[m_nyPP]<m_numLines[m_nyPP]-1; ++posL[m_nyPP]){
            posR[m_nyP]  = posL[m_nyP];
            posR[m_nyPP] = posL[m_nyPP];
            tmp_Iim[0] = curr_im[0][posL[0]][posL[1]][posL[2]];
            tmp_Iim[1] = curr_im[1][posL[0]][posL[1]][posL[2]];
            tmp_Iim[2] = curr_im[2][posL[0]][posL[1]][posL[2]];
            tmp_Ire[0] = m_Eng->curr[0][posL[0]][posL[1]][posL[2]];
            tmp_Ire[1] = m_Eng->curr[1][posL[0]][posL[1]][posL[2]];
            tmp_Ire[2] = m_Eng->curr[2][posL[0]][posL[1]][posL[2]];
            tmp_Uim[0] = volt_im[0][posL[0]][posL[1]][posL[2]];
            tmp_Uim[1] = volt_im[1][posL[0]][posL[1]][posL[2]];
            tmp_Uim[2] = volt_im[2][posL[0]][posL[1]][posL[2]];
            tmp_Ure[0] = m_Eng->volt[0][posL[0]][posL[1]][posL[2]];
            tmp_Ure[1] = m_Eng->volt[1][posL[0]][posL[1]][posL[2]];
            tmp_Ure[2] = m_Eng->volt[2][posL[0]][posL[1]][posL[2]];
            // apply phases to imaginary and real parts of the voltage/current on the left pbc border
            curr_im[0][posL[0]][posL[1]][posL[2]] = cosinus*curr_im[0][posR[0]][posR[1]][posR[2]]+sinus*m_Eng->GetCurr(0, posR);
            curr_im[1][posL[0]][posL[1]][posL[2]] = cosinus*curr_im[1][posR[0]][posR[1]][posR[2]]+sinus*m_Eng->GetCurr(1, posR);
            curr_im[2][posL[0]][posL[1]][posL[2]] = cosinus*curr_im[2][posR[0]][posR[1]][posR[2]]+sinus*m_Eng->GetCurr(2, posR);
            volt_im[0][posL[0]][posL[1]][posL[2]] = cosinus*volt_im[0][posR[0]][posR[1]][posR[2]]+sinus*m_Eng->GetVolt(0, posR);
            volt_im[1][posL[0]][posL[1]][posL[2]] = cosinus*volt_im[1][posR[0]][posR[1]][posR[2]]+sinus*m_Eng->GetVolt(1, posR);
            volt_im[2][posL[0]][posL[1]][posL[2]] = cosinus*volt_im[2][posR[0]][posR[1]][posR[2]]+sinus*m_Eng->GetVolt(2, posR);

            m_Eng->SetCurr(0, posL, cosinus*m_Eng->GetCurr(0,posR) - sinus*curr_im[0][posR[0]][posR[1]][posR[2]]);
            m_Eng->SetCurr(1, posL, cosinus*m_Eng->GetCurr(1,posR) - sinus*curr_im[1][posR[0]][posR[1]][posR[2]]);
            m_Eng->SetCurr(2, posL, cosinus*m_Eng->GetCurr(2,posR) - sinus*curr_im[2][posR[0]][posR[1]][posR[2]]);
            m_Eng->SetVolt(0, posL, cosinus*m_Eng->GetVolt(0,posR) - sinus*volt_im[0][posR[0]][posR[1]][posR[2]]);
            m_Eng->SetVolt(1, posL, cosinus*m_Eng->GetVolt(1,posR) - sinus*volt_im[1][posR[0]][posR[1]][posR[2]]);
            m_Eng->SetVolt(2, posL, cosinus*m_Eng->GetVolt(2,posR) - sinus*volt_im[2][posR[0]][posR[1]][posR[2]]);
            // apply phases to imaginary and real parts of the voltage/current on the right pbc border
            curr_im[0][posR[0]][posR[1]][posR[2]] = cosinus*tmp_Iim[0]-sinus*tmp_Ire[0];
            curr_im[1][posR[0]][posR[1]][posR[2]] = cosinus*tmp_Iim[1]-sinus*tmp_Ire[1];
            curr_im[2][posR[0]][posR[1]][posR[2]] = cosinus*tmp_Iim[2]-sinus*tmp_Ire[2];
            volt_im[0][posR[0]][posR[1]][posR[2]] = cosinus*tmp_Uim[0]-sinus*tmp_Ure[0];
            volt_im[1][posR[0]][posR[1]][posR[2]] = cosinus*tmp_Uim[1]-sinus*tmp_Ure[1];
            volt_im[2][posR[0]][posR[1]][posR[2]] = cosinus*tmp_Uim[2]-sinus*tmp_Ure[2];

            m_Eng->SetCurr(0, posR, cosinus*tmp_Ire[0] + sinus*tmp_Iim[0]);
            m_Eng->SetCurr(1, posR, cosinus*tmp_Ire[1] + sinus*tmp_Iim[1]);
            m_Eng->SetCurr(2, posR, cosinus*tmp_Ire[2] + sinus*tmp_Iim[2]);
            m_Eng->SetVolt(0, posR, cosinus*tmp_Ure[0] + sinus*tmp_Uim[0]);
            m_Eng->SetVolt(1, posR, cosinus*tmp_Ure[1] + sinus*tmp_Uim[1]);
            m_Eng->SetVolt(2, posR, cosinus*tmp_Ure[2] + sinus*tmp_Uim[2]);
        }
    }

};

void Engine_Ext_Pbc::DoPreVoltageUpdates(int threadID){
    for(unsigned int dir = 0; dir<3; ++dir){
        if (direction_is_pbc[dir])
            Apply_Phases_to_dir(dir);
    }
};
void Engine_Ext_Pbc::Apply2Voltages(int threadID){};

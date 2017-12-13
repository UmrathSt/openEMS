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
    k_pbc = m_Op_Pbc->k_pbc;
    sin_kx = sin(k_pbc[0]);
    cos_kx = cos(k_pbc[0]);
    sin_ky = sin(k_pbc[1]);
    cos_ky = cos(k_pbc[1]);
    cos_kxy = cos(k_pbc[1] + k_pbc[0]);
    sin_kxy = sin(k_pbc[1] + k_pbc[0]);
    cout << "engine_ext_pbc.cpp: called the constructor" << endl;
    cout << "engine_ext_pbc.cpp: kx,ky,kz=" << k_pbc[0] << "," << k_pbc[1] << "," << k_pbc[2] << endl;
    m_numLines[0] = m_Op_Pbc->m_numLines[0];
    m_numLines[1] = m_Op_Pbc->m_numLines[1];
    m_numLines[2] = m_Op_Pbc->m_numLines[2];
    maxX = m_numLines[0]-1;
    maxY = m_numLines[1]-1;
    volt_im = Create_N_3DArray<FDTD_FLOAT>(m_numLines); // imaginary parts of
    curr_im = Create_N_3DArray<FDTD_FLOAT>(m_numLines); // voltage/current as (3, Nx, Ny, Nz) array of floats
    SetNumberOfThreads(1);

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

void Engine_Ext_Pbc::Apply2Voltages(){

    // use peridic boundary condition to obtain the current values which are outside the FDTD domain
    pos[0]=0;
    for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2]){// moved the case y = pos[1] = 0 outside the loop
        pos[1] = 0;
        curr_im_outside[0] = cos_kxy*curr_im[0][maxX][maxY][pos[2]] + sin_kxy*m_Eng->GetCurr(0,maxX,maxY,pos[2]);
        curr_im_outside[1] = cos_kxy*curr_im[1][maxX][maxY][pos[2]] + sin_kxy*m_Eng->GetCurr(1,maxX,maxY,pos[2]);
        curr_im_outside[2] = cos_kxy*curr_im[2][maxX][maxY][pos[2]] + sin_kxy*m_Eng->GetCurr(2,maxX,maxY,pos[2]);
        curr_outside[0] = cos_kxy*m_Eng->GetCurr(0,maxX,maxY,pos[2]) - sin_kxy*curr_im[0][maxX][maxY][pos[2]];
        curr_outside[1] = cos_kxy*m_Eng->GetCurr(1,maxX,maxY,pos[2]) - sin_kxy*curr_im[1][maxX][maxY][pos[2]];
        curr_outside[2] = cos_kxy*m_Eng->GetCurr(2,maxX,maxY,pos[2]) - sin_kxy*curr_im[2][maxX][maxY][pos[2]];
        volt_im[0][maxX][maxY][pos[2]] *= m_Op_Pbc->GetVV(0,maxX,maxY,pos[2]);
        volt_im[0][maxX][maxY][pos[2]] += m_Op_Pbc->GetVI(0,maxX,maxY,pos[2]) * (curr_im[2][maxX][maxY][pos[2]] - curr_im_outside[2] - curr_im[1][maxX][maxY][pos[2]] + curr_im_outside[1]);
        m_Eng->SetVolt(0,maxX,maxY,pos[2], m_Eng->GetVolt(0,maxX,maxY,pos[2]) +m_Op_Pbc->GetVI(0,maxX,maxY,pos[2]) * (
                    m_Eng->GetCurr(2,maxX,maxY,pos[2]) - curr_outside[2] - m_Eng->GetCurr(1,maxX,maxY,pos[2]) + curr_outside[1]));
        //for y
        volt_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(1,pos[0],pos[1],pos[2]);
        volt_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(1,pos[0],pos[1],pos[2]) * (curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im_outside[0] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im_outside[2]);
        m_Eng->SetVolt(1,maxX,maxY,pos[2], m_Eng->GetVolt(1,maxX,maxY,pos[2]) +m_Op_Pbc->GetVI(1,maxX,maxY,pos[2]) * (
                    m_Eng->GetCurr(0,maxX,maxY,pos[2]) - curr_outside[0] - m_Eng->GetCurr(2,maxX,maxY,pos[2]) + curr_outside[2]));
        //for z
        volt_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(2,pos[0],pos[1],pos[2]);
        volt_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(2,pos[0],pos[1],pos[2]) * (curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im_outside[1] - curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im_outside[0]);
        m_Eng->SetVolt(2,maxX,maxY,pos[2], m_Eng->GetVolt(2,maxX,maxY,pos[2]) +m_Op_Pbc->GetVI(2,maxX,maxY,pos[2]) * (
                    m_Eng->GetCurr(1,maxX,maxY,pos[2]) - curr_outside[1] - m_Eng->GetCurr(0,maxX,maxY,pos[2]) + curr_outside[0]));
        for (pos[1]=1; pos[1]<m_numLines[1]; ++pos[1]){
            curr_im_outside[0] = cos_kx*volt_im[0][maxX][pos[1]][pos[2]] + sin_kx*m_Eng->GetVolt(0,maxX,pos[1],pos[2]);
            curr_im_outside[1] = cos_kx*volt_im[1][maxX][pos[1]][pos[2]] + sin_kx*m_Eng->GetVolt(1,maxX,pos[1],pos[2]);
            curr_im_outside[2] = cos_kx*volt_im[2][maxX][pos[1]][pos[2]] + sin_kx*m_Eng->GetVolt(2,maxX,pos[1],pos[2]);
            curr_outside[0] = cos_kx*m_Eng->GetVolt(0,maxX,pos[1],pos[2]) - sin_kx*volt_im[0][maxX][pos[1]][pos[2]];
            curr_outside[1] = cos_kx*m_Eng->GetVolt(1,maxX,pos[1],pos[2]) - sin_kx*volt_im[1][maxX][pos[1]][pos[2]];
            curr_outside[2] = cos_kx*m_Eng->GetVolt(2,maxX,pos[1],pos[2]) - sin_kx*volt_im[2][maxX][pos[1]][pos[2]];
            //do the updates here
            //for x
            volt_im[0][maxX][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(0,maxX,pos[1],pos[2]);
            volt_im[0][maxX][pos[1]][pos[2]] += m_Op_Pbc->GetVI(0,maxX,pos[1],pos[2]) * (curr_im[2][maxX][pos[1]][pos[2]] - curr_im_outside[2] - curr_im[1][maxX][pos[1]][pos[2]] + curr_im_outside[1]);
            m_Eng->SetVolt(0,maxX,pos[1],pos[2], m_Eng->GetVolt(0,maxX,pos[1],pos[2]) +m_Op_Pbc->GetVI(0,maxX,pos[1],pos[2]) * (
                        m_Eng->GetCurr(2,maxX,pos[1],pos[2]) - curr_outside[2] - m_Eng->GetCurr(1,maxX,pos[1],pos[2]) + curr_outside[1]));
            //for y
            volt_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(1,pos[0],pos[1],pos[2]);
            volt_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(1,pos[0],pos[1],pos[2]) * (curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im_outside[0] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im_outside[2]);
            m_Eng->SetVolt(1,maxX,pos[1],pos[2], m_Eng->GetVolt(1,maxX,pos[1],pos[2]) +m_Op_Pbc->GetVI(1,maxX,pos[1],pos[2]) * (
                        m_Eng->GetCurr(0,maxX,pos[1],pos[2]) - curr_outside[0] - m_Eng->GetCurr(2,maxX,pos[1],pos[2]) + curr_outside[2]));
            //for z
            volt_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(2,pos[0],pos[1],pos[2]);
            volt_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(2,pos[0],pos[1],pos[2]) * (curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im_outside[1] - curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im_outside[0]);
            m_Eng->SetVolt(2,maxX,pos[1],pos[2], m_Eng->GetVolt(2,maxX,pos[1],pos[2]) +m_Op_Pbc->GetVI(2,maxX,pos[1],pos[2]) * (
                        m_Eng->GetCurr(1,maxX,pos[1],pos[2]) - curr_outside[1] - m_Eng->GetCurr(0,maxX,pos[1],pos[2]) + curr_outside[0]));
        }

        }

    //
    for (pos[0]=1; pos[0]<m_numLines[0]; ++pos[0])
    {
        shift[0]=pos[0];
        for (pos[1]=1; pos[1]<m_numLines[1]; ++pos[1])
        {
            shift[1]=pos[1]; // 1 for all but the zeroth gird line, since shift is of type "bool"
            for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
            {
                shift[2]=pos[2]; // 1 for all but the zeroth gird line
                //do the updates here
                //for x
                volt_im[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(0,pos[0],pos[1],pos[2]);
                volt_im[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(0,pos[0],pos[1],pos[2]) * (curr_im[2][pos[0]][pos[1]][pos[2]] - curr_im[2][pos[0]][pos[1]-shift[1]][pos[2]] - curr_im[1][pos[0]][pos[1]][pos[2]] + curr_im[1][pos[0]][pos[1]][pos[2]-shift[2]]);

                //for y
                volt_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(1,pos[0],pos[1],pos[2]);
                volt_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(1,pos[0],pos[1],pos[2]) * (curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im[0][pos[0]][pos[1]][pos[2]-shift[2]] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im[2][pos[0]-shift[0]][pos[1]][pos[2]]);

                //for z
                volt_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVV(2,pos[0],pos[1],pos[2]);
                volt_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVI(2,pos[0],pos[1],pos[2]) * (curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im[1][pos[0]-shift[0]][pos[1]][pos[2]]- curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im[0][pos[0]][pos[1]-shift[1]][pos[2]]);
            }
        }
        ++pos[0];
    }

}


void Engine_Ext_Pbc::Apply2Current(){
    // to update the currents at the lower X and Y boundary we need current values outside the domain
    // Those are available from currents at the opposite side
    pos[0] = maxX;
    for (pos[2]=0; pos[2]<m_numLines[2]-1; ++pos[2]){
        pos[1] = maxY;
        volt_im_outside[0] = cos_kxy*volt_im[0][0][0][pos[2]] - sin_kxy*m_Eng->GetVolt(0,0,0,pos[2]);
        volt_im_outside[1] = cos_kxy*volt_im[1][0][0][pos[2]] - sin_kxy*m_Eng->GetVolt(1,0,0,pos[2]);
        volt_im_outside[2] = cos_kxy*volt_im[2][0][0][pos[2]] - sin_kxy*m_Eng->GetVolt(2,0,0,pos[2]);
        volt_outside[0] = cos_kxy*m_Eng->GetVolt(0,0,0,pos[2]) + sin_kxy*volt_im[0][0][0][pos[2]];
        volt_outside[1] = cos_kxy*m_Eng->GetVolt(1,0,0,pos[2]) + sin_kxy*volt_im[1][0][0][pos[2]];
        volt_outside[2] = cos_kxy*m_Eng->GetVolt(2,0,0,pos[2]) + sin_kxy*volt_im[2][0][0][pos[2]];
        curr_im[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(0,pos[0],pos[1],pos[2]);
        curr_im[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(0,pos[0],pos[1],pos[2]) * (volt_im[2][pos[0]][pos[1]][pos[2]] - volt_im_outside[2] - curr_im[1][pos[0]][pos[1]][pos[2]] + curr_im_outside[1]);
        m_Eng->SetCurr(0,pos[0],pos[1],pos[2], m_Eng->GetVolt(0,0,0,pos[2]) +m_Op_Pbc->GetII(0,0,0,pos[2]) * (
                    m_Eng->GetVolt(2,pos[0],pos[1],pos[2]) - volt_outside[2] - m_Eng->GetVolt(1,pos[0],pos[1],pos[2]) + volt_outside[1]));
        //for y
        curr_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(1,pos[0],pos[1],pos[2]);
        curr_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(1,pos[0],pos[1],pos[2]) * (curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im_outside[0] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im_outside[2]);
        m_Eng->SetCurr(1,pos[0],pos[1],pos[2], m_Eng->GetVolt(1,0,0,pos[2]) +m_Op_Pbc->GetIV(1,0,0,pos[2]) * (
                    m_Eng->GetVolt(0,pos[0],pos[1],pos[2]) - volt_outside[0] - m_Eng->GetVolt(2,pos[0],pos[1],pos[2]) + volt_outside[2]));
        //for z
        curr_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(2,pos[0],pos[1],pos[2]);
        curr_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(2,pos[0],pos[1],pos[2]) * (curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im_outside[1] - curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im_outside[0]);
        m_Eng->SetCurr(2,pos[0],pos[1],pos[2], m_Eng->GetVolt(2,0,0,pos[2]) +m_Op_Pbc->GetIV(2,0,0,pos[2]) * (
                    m_Eng->GetVolt(1,pos[0],pos[1],pos[2]) - volt_outside[1] - m_Eng->GetVolt(0,pos[0],pos[1],pos[2]) + volt_outside[0]));
        for (pos[1]=0; pos[1]<m_numLines[1]-1; ++pos[1]){
            volt_im_outside[0] = cos_kx*volt_im[0][0][pos[1]][pos[2]] + sin_kx*m_Eng->GetVolt(0,0,pos[1],pos[2]);
            volt_im_outside[1] = cos_kx*volt_im[1][0][pos[1]][pos[2]] + sin_kx*m_Eng->GetVolt(1,0,pos[1],pos[2]);
            volt_im_outside[2] = cos_kx*volt_im[2][0][pos[1]][pos[2]] + sin_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]);
            volt_outside[0] = cos_kx*m_Eng->GetVolt(0,0,pos[1],pos[2]) - sin_kx*volt_im[0][0][pos[1]][pos[2]];
            volt_outside[1] = cos_kx*m_Eng->GetVolt(1,0,pos[1],pos[2]) - sin_kx*volt_im[1][0][pos[1]][pos[2]];
            volt_outside[2] = cos_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]) - sin_kx*volt_im[2][0][pos[1]][pos[2]];
            //do the updates here
            //for x
            curr_im[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(0,pos[0],pos[1],pos[2]);

            curr_im[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(0,pos[0],pos[1],pos[2]) * (volt_im[2][0][pos[1]][pos[2]] - volt_im_outside[2] - volt_im[1][pos[0]][pos[1]][pos[2]] + volt_im_outside[1]);
            m_Eng->SetCurr(0,pos[0],pos[1],pos[2], m_Eng->GetCurr(0,pos[0],pos[1],pos[2]) +m_Op_Pbc->GetIV(0,0,pos[1],pos[2]) * (
                        m_Eng->GetVolt(2,pos[0],pos[1],pos[2]) - volt_outside[2] - m_Eng->GetVolt(1,pos[0],pos[1],pos[2]) + volt_outside[1]));
            //for y
            curr_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(1,pos[0],pos[1],pos[2]);
            curr_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(1,pos[0],pos[1],pos[2]) * (volt_im[0][pos[0]][pos[1]][pos[2]] - volt_im_outside[0] - volt_im[2][pos[0]][pos[1]][pos[2]] + volt_im_outside[2]);
            m_Eng->SetCurr(1,pos[0],pos[1],pos[2], m_Eng->GetCurr(1,pos[0],pos[1],pos[2]) +m_Op_Pbc->GetIV(1,pos[0],pos[1],pos[2]) * (
                        m_Eng->GetVolt(0,pos[0],pos[1],pos[2]) - volt_outside[0] - m_Eng->GetVolt(2,pos[0],pos[1],pos[2]) + volt_outside[2]));
            //for z
            curr_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(2,pos[0],pos[1],pos[2]);
            curr_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(2,pos[0],pos[1],pos[2]) * (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im_outside[1] - volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im_outside[0]);
            m_Eng->SetCurr(2,pos[0],pos[1],pos[2], m_Eng->GetCurr(2,pos[0],pos[1],pos[2]) +m_Op_Pbc->GetIV(2,pos[0],pos[1],pos[2]) * (
                        m_Eng->GetVolt(1,pos[0],pos[1],pos[2]) - volt_outside[1] - m_Eng->GetVolt(0,pos[0],pos[1],pos[2]) + volt_outside[0]));
        }}

    for (pos[0]=0; pos[0]<m_numLines[0]-1; ++pos[0])
    {
        shift[0] = pos[0];
        for (pos[1]=0; pos[1]<m_numLines[1]-1; ++pos[1])
        {
            shift[1] = pos[1];
            for (pos[2]=0; pos[2]<m_numLines[2]-1; ++pos[2])
            {
                shift[2] = pos[2];
                //do the updates here
                //for x
                curr_im[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(0,pos[0],pos[1],pos[2]);
                curr_im[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(0,pos[0],pos[1],pos[2]) * (volt_im[2][pos[0]][pos[1]][pos[2]] - volt_im[2][pos[0]][pos[1]+1][pos[2]] - volt_im[1][pos[0]][pos[1]][pos[2]] + volt_im[1][pos[0]][pos[1]][pos[2]+1]);

                //for y
                curr_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(1,pos[0],pos[1],pos[2]);
                curr_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(1,pos[0],pos[1],pos[2]) * (volt_im[0][pos[0]][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]+1] - volt_im[2][pos[0]][pos[1]][pos[2]] + volt_im[2][pos[0]+1][pos[1]][pos[2]]);
                //for z
                curr_im[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetII(2,pos[0],pos[1],pos[2]);
                curr_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(2,pos[0],pos[1],pos[2]) * (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im[1][pos[0]+1][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im[0][pos[0]][pos[1]+1][pos[2]]);

            }
        }
        ++pos[0];
    }
};




void Engine_Ext_Pbc::DoPreVoltageUpdates(int threadID)
{
};


void Engine_Ext_Pbc::DoPostCurrentUpdates(int threadID){
    //soft current excitation here (E-field excite)
    int exc_pos;
    unsigned int ny;
    unsigned int pos[3];
    int numTS = m_Eng->GetNumberOfTimesteps();
    unsigned int length = m_Op_Pbc->m_Exc->GetLength();
    FDTD_FLOAT* exc_curr_sin =  m_Op_Pbc->m_Exc->GetCurrentSignal_s();
    FDTD_FLOAT* exc_curr_cos =  m_Op_Pbc->m_Exc->GetCurrentSignal(); // the signals time-dependence
    int p = numTS+1;
    if (m_Op_Pbc->m_Exc->GetSignalPeriod()>0)
        p = int(m_Op_Pbc->m_Exc->GetSignalPeriod()/m_Op_Pbc->m_Exc->GetTimestep());
    //switch for different engine types to access faster inline engine functions
    switch (m_Eng->GetType())
    {
    case Engine::BASIC:
        {
            for (unsigned int n=0; n<m_Op_Pbc->Curr_Count; ++n)
            {
                exc_pos = numTS - (int)m_Op_Pbc->Curr_delay[n];
                exc_pos *= (exc_pos>0);
                exc_pos %= p;
                exc_pos *= (exc_pos<(int)length);
                ny = m_Op_Pbc->Curr_dir[n];
                pos[0]=m_Op_Pbc->Curr_index[0][n];
                pos[1]=m_Op_Pbc->Curr_index[1][n];
                pos[2]=m_Op_Pbc->Curr_index[2][n];
                m_Eng->Engine::SetCurr(ny,pos, m_Eng->Engine::GetCurr(ny,pos)  + m_Op_Pbc->Curr_amp_sin[n]*exc_curr_sin[exc_pos]
                                                                               + m_Op_Pbc->Curr_amp_cos[n]*exc_curr_cos[exc_pos]);
            }
            break;
        }
     default:
    {
        break;
    }
    }
};
void Engine_Ext_Pbc::DoPostVoltageUpdates(int threadID){
    //soft voltage excitation here (E-field excite)
    int exc_pos;
    unsigned int ny;
    unsigned int pos[3];
    int numTS = m_Eng->GetNumberOfTimesteps();
    unsigned int length = m_Op_Pbc->m_Exc->GetLength();
    FDTD_FLOAT* exc_volt_sin =  m_Op_Pbc->m_Exc->GetVoltageSignal_s();
    FDTD_FLOAT* exc_volt_cos =  m_Op_Pbc->m_Exc->GetVoltageSignal();

    int p = numTS+1;
    if (m_Op_Pbc->m_Exc->GetSignalPeriod()>0)
        p = int(m_Op_Pbc->m_Exc->GetSignalPeriod()/m_Op_Pbc->m_Exc->GetTimestep());

    //switch for different engine types to access faster inline engine functions
    switch (m_Eng->GetType())
    {
    case Engine::BASIC:
        {
            for (unsigned int n=0; n<m_Op_Pbc->Volt_Count; ++n)
            {
                exc_pos = numTS - (int)m_Op_Pbc->Volt_delay[n];
                exc_pos *= (exc_pos>0);
                exc_pos %= p;
                exc_pos *= (exc_pos<(int)length);
                ny = m_Op_Pbc->Volt_dir[n];
                pos[0]=m_Op_Pbc->Volt_index[0][n];
                pos[1]=m_Op_Pbc->Volt_index[1][n];
                pos[2]=m_Op_Pbc->Volt_index[2][n];
                m_Eng->Engine::SetVolt(ny,pos, m_Eng->Engine::GetVolt(ny,pos) + m_Op_Pbc->Volt_amp_sin[n]*exc_volt_sin[exc_pos]
                                                                              + m_Op_Pbc->Volt_amp_cos[n]*exc_volt_cos[exc_pos]);
            }
            break;
        }
     default:
    {
        break;
    }
    }
};

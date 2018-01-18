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
    if(!(m_Op_Pbc->dir_is_pbc[0] && m_Op_Pbc->dir_is_pbc[1] && m_Op_Pbc->dir_is_pbc[2] && m_Op_Pbc->dir_is_pbc[3])){
        cerr << "engine_ext_pbc.cpp: Error +-X and +-Y directions are expected to be PBC." << endl;
        exit(-3);
    }
    pbc_phase = m_Op_Pbc->pbc_phase;
    sin_kx = sin(pbc_phase[0]);
    cos_kx = cos(pbc_phase[0]);
    sin_ky = sin(pbc_phase[1]);
    cos_ky = cos(pbc_phase[1]);
    std::cout << "engine_ext_pbc.cpp: called the constructor" << endl;
    std::cout << "engine_ext_pbc.cpp: kx,ky,kz=" << pbc_phase[0] << "," << pbc_phase[1] << "," << pbc_phase[2] << std::endl;
    m_numLines[0] = m_Op_Pbc->m_numLines[0];
    m_numLines[1] = m_Op_Pbc->m_numLines[1];
    m_numLines[2] = m_Op_Pbc->m_numLines[2];
    volt_im = Create_N_3DArray<FDTD_FLOAT>(m_numLines); // imaginary parts of
    curr_im = Create_N_3DArray<FDTD_FLOAT>(m_numLines); // voltage/current as (3, Nx, Ny, Nz) array of floats
    maxX = m_numLines[0]-1;
    maxY = m_numLines[1]-1;

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



void Engine_Ext_Pbc::DoCurrPhaseUpdates(){
    bool shift;

    for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
    {
        pos[0]=maxX;
        // moved the case x = pos[maxX] = maxX outside the loop
        // now treat the YZ-plane at the minimum x-coordinate and get the missing currents from pos[0] = maxX
        for (pos[1]=0; pos[1]<m_numLines[1]-1; ++pos[1]){
            volt_im_outside[1] = cos_kx*volt_im[1][0][pos[1]][pos[2]]     + sin_kx*m_Eng->GetVolt(1,0,pos[1],pos[2]);
            volt_im_outside[2] = cos_kx*volt_im[2][0][pos[1]][pos[2]]     + sin_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]);
            volt_outside[1] =    cos_kx*m_Eng->GetVolt(1,0,pos[1],pos[2]) - sin_kx*volt_im[1][0][pos[1]][pos[2]];
            volt_outside[2] =    cos_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]) - sin_kx*volt_im[2][0][pos[1]][pos[2]];

            m_Eng->curr[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(0,pos[0],pos[1],pos[2]);
            m_Eng->curr[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(0,pos[0],pos[1],pos[2])* ( m_Eng->volt[2][pos[0]][pos[1]][pos[2]] - m_Eng->volt[2][pos[0]][pos[1]+1][pos[2]] - m_Eng->volt[1][pos[0]][pos[1]][pos[2]] + m_Eng->volt[1][pos[0]][pos[1]][pos[2]+1]);
            curr_im[0][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(0,pos[0],pos[1],pos[2]);
            curr_im[0][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(0,pos[0],pos[1],pos[2]) * (volt_im[2][pos[0]][pos[1]][pos[2]] - volt_im[2][pos[0]][pos[1]+1][pos[2]] - volt_im[1][pos[0]][pos[1]][pos[2]] + volt_im[1][pos[0]][pos[1]][pos[2]+1]);

            //for y
            m_Eng->curr[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(1,pos[0],pos[1],pos[2]);
            m_Eng->curr[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(1,pos[0],pos[1],pos[2])* ( m_Eng->volt[0][pos[0]][pos[1]][pos[2]] - m_Eng->volt[0][pos[0]][pos[1]][pos[2]+1] - m_Eng->volt[2][pos[0]][pos[1]][pos[2]] + volt_outside[2]);
            curr_im[1][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(1,pos[0],pos[1],pos[2]);
            curr_im[1][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(1,pos[0],pos[1],pos[2]) * (volt_im[0][pos[0]][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]+1] - volt_im[2][pos[0]][pos[1]][pos[2]] + volt_im_outside[2]);

            //for z
            m_Eng->curr[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(2,pos[0],pos[1],pos[2]);
            m_Eng->curr[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(2,pos[0],pos[1],pos[2])* ( m_Eng->volt[1][pos[0]][pos[1]][pos[2]] - volt_outside[1]- m_Eng->volt[0][pos[0]][pos[1]][pos[2]] + m_Eng->volt[0][pos[0]][pos[1]+1][pos[2]]);
            curr_im[2][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(2,pos[0],pos[1],pos[2]);
            curr_im[2][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(2,pos[0],pos[1],pos[2]) * (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im_outside[1]- volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im[0][pos[0]][pos[1]+1][pos[2]]);

        }
        pos[1] = maxY;
        // moved the case y = pos[1] = maxY outside the loop
        // now treat the XZ-plane at the minimum x-coordinate and get the missing currents from pos[0] = 0
        for (pos[0]=0; pos[0]<m_numLines[0]-1; ++pos[0]){
            volt_im_outside[0] = cos_ky*volt_im[0][pos[0]][0][pos[2]]     + sin_ky*m_Eng->GetVolt(0,pos[0],0,pos[2]);
            volt_im_outside[2] = cos_ky*volt_im[2][pos[0]][0][pos[2]]     + sin_ky*m_Eng->GetVolt(2,pos[0],0,pos[2]);
            volt_outside[0] =    cos_ky*m_Eng->GetVolt(0,pos[0],0,pos[2]) - sin_ky*volt_im[0][pos[0]][0][pos[2]];
            volt_outside[2] =    cos_ky*m_Eng->GetVolt(2,pos[0],0,pos[2]) - sin_ky*volt_im[2][pos[0]][0][pos[2]];

            m_Eng->curr[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(0,pos[0],pos[1],pos[2]);
            m_Eng->curr[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(0,pos[0],pos[1],pos[2])* ( m_Eng->volt[2][pos[0]][pos[1]][pos[2]] - volt_outside[2]- m_Eng->volt[1][pos[0]][pos[1]][pos[2]] + m_Eng->volt[1][pos[0]][pos[1]][pos[2]+1]);
            curr_im[0][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(0,pos[0],pos[1],pos[2]);
            curr_im[0][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(0,pos[0],pos[1],pos[2]) * (volt_im[2][pos[0]][pos[1]][pos[2]] - volt_im_outside[2] - volt_im[1][pos[0]][pos[1]][pos[2]] + volt_im[1][pos[0]][pos[1]][pos[2]+1]);
            //for y
            m_Eng->curr[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(1,pos[0],pos[1],pos[2]);
            m_Eng->curr[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(1,pos[0],pos[1],pos[2])* ( m_Eng->volt[0][pos[0]][pos[1]][pos[2]] - m_Eng->volt[0][pos[0]][pos[1]][pos[2]+1] - m_Eng->volt[2][pos[0]][pos[1]][pos[2]] + m_Eng->volt[2][pos[0]+1][pos[1]][pos[2]]);
            curr_im[1][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(1,pos[0],pos[1],pos[2]);
            curr_im[1][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(1,pos[0],pos[1],pos[2]) * (volt_im[0][pos[0]][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]+1] - volt_im[2][pos[0]][pos[1]][pos[2]] + volt_im[2][pos[0]+1][pos[1]][pos[2]]);
            //for z
            m_Eng->curr[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(2,pos[0],pos[1],pos[2]);
            m_Eng->curr[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(2,pos[0],pos[1],pos[2])* ( m_Eng->volt[1][pos[0]][pos[1]][pos[2]] - m_Eng->volt[1][pos[0]+1][pos[1]][pos[2]] - m_Eng->volt[0][pos[0]][pos[1]][pos[2]] + volt_outside[0]);
            curr_im[2][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(2,pos[0],pos[1],pos[2]);
            curr_im[2][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(2,pos[0],pos[1],pos[2]) * (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im[1][pos[0]+1][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im_outside[0]);
                    }
        pos[0] = maxX;
        pos[1] = maxY;
        volt_outside[2] =    cos_ky*m_Eng->GetVolt(2,pos[0],0,pos[2]) - sin_ky*volt_im[2][pos[0]][0][pos[2]];
        volt_im_outside[2] = cos_ky*volt_im[2][pos[0]][0][pos[2]]     + sin_ky*m_Eng->GetVolt(2,pos[0],0,pos[2]);
        // for x updates
        m_Eng->curr[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(0,pos[0],pos[1],pos[2]);
        m_Eng->curr[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(0,pos[0],pos[1],pos[2])* ( m_Eng->volt[2][pos[0]][pos[1]][pos[2]] - volt_outside[2]- m_Eng->volt[1][pos[0]][pos[1]][pos[2]] + m_Eng->volt[1][pos[0]][pos[1]][pos[2]+1]);
        curr_im[0][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(0,pos[0],pos[1],pos[2]);
        curr_im[0][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(0,pos[0],pos[1],pos[2])* (volt_im[2][pos[0]][pos[1]][pos[2]] - volt_im_outside[2] - volt_im[1][pos[0]][pos[1]][pos[2]] + volt_im[1][pos[0]][pos[1]][pos[2]+1]);
        volt_outside[2] =    cos_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]) - sin_kx*volt_im[2][0][pos[1]][pos[2]];
        volt_im_outside[2] = cos_kx*volt_im[2][0][pos[1]][pos[2]]  + sin_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]);
        // for y updates
        m_Eng->curr[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(1,pos[0],pos[1],pos[2]);
        m_Eng->curr[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(1,pos[0],pos[1],pos[2])* ( m_Eng->volt[0][pos[0]][pos[1]][pos[2]] - m_Eng->volt[0][pos[0]][pos[1]][pos[2]+1] - m_Eng->volt[2][pos[0]][pos[1]][pos[2]] + volt_outside[2]);
        curr_im[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(1,pos[0],pos[1],pos[2]);
        curr_im[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(1,pos[0],pos[1],pos[2])* (volt_im[0][pos[0]][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]+1] - volt_im[2][pos[0]][pos[1]][pos[2]] + volt_im_outside[2]);
        volt_outside[0] =    cos_ky*m_Eng->GetVolt(0,pos[0],0,pos[2]) - sin_ky*volt_im[0][pos[0]][0][pos[2]];
        volt_outside[1] =    cos_kx*m_Eng->GetVolt(1,0,pos[1],pos[2]) - sin_kx*volt_im[1][0][pos[1]][pos[2]];
        volt_im_outside[0] = cos_ky*volt_im[2][pos[0]][0][pos[2]]     + sin_ky*m_Eng->GetVolt(2,pos[0],0,pos[2]);
        volt_im_outside[1] = cos_kx*volt_im[2][0][pos[1]][pos[2]]     + sin_kx*m_Eng->GetVolt(2,0,pos[1],pos[2]);
        //for z updates
        m_Eng->curr[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetIIedge(2,pos[0],pos[1],pos[2]);
        m_Eng->curr[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIVedge(2,pos[0],pos[1],pos[2])* ( m_Eng->volt[1][pos[0]][pos[1]][pos[2]] - volt_outside[1] - m_Eng->volt[0][pos[0]][pos[1]][pos[2]] + volt_outside[0]);
        curr_im[2][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetIIedge(2,pos[0],pos[1],pos[2]);
        curr_im[2][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetIVedge(2,pos[0],pos[1],pos[2])* (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im_outside[1] - volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im_outside[0]);
    }

};




void Engine_Ext_Pbc::DoVoltPhaseUpdates()
{
    bool shift;

    for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
    {
        pos[0]=0;
        shift = pos[2];
        // moved the case x = pos[0] = 0 outside the loop
        // now treat the YZ-plane at the minimum x-coordinate and get the missing currents from pos[0] = maxX
        for (pos[1]=1; pos[1]<m_numLines[1]; ++pos[1]){
            curr_im_outside[1] = cos_kx*curr_im[1][maxX][pos[1]][pos[2]]     - sin_kx*m_Eng->GetCurr(1,maxX,pos[1],pos[2]);
            curr_im_outside[2] = cos_kx*curr_im[2][maxX][pos[1]][pos[2]]     - sin_kx*m_Eng->GetCurr(2,maxX,pos[1],pos[2]);
            curr_outside[1] =    cos_kx*m_Eng->curr[1][maxX][pos[1]][pos[2]] + sin_kx*curr_im[1][maxX][pos[1]][pos[2]];
            curr_outside[2] =    cos_kx*m_Eng->curr[2][maxX][pos[1]][pos[2]] + sin_kx*curr_im[2][maxX][pos[1]][pos[2]];

            m_Eng->volt[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(0,pos[0],pos[1],pos[2]);
            m_Eng->volt[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(0,pos[0],pos[1],pos[2])* ( m_Eng->curr[2][pos[0]][pos[1]][pos[2]] - m_Eng->curr[2][pos[0]][pos[1]-1][pos[2]] - m_Eng->curr[1][pos[0]][pos[1]][pos[2]] + m_Eng->curr[1][pos[0]][pos[1]][pos[2]-shift]);
            volt_im[0][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(0,pos[0],pos[1],pos[2]);
            volt_im[0][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(0,pos[0],pos[1],pos[2])* ( curr_im[2][pos[0]][pos[1]][pos[2]] - curr_im[2][pos[0]][pos[1]-1][pos[2]] - curr_im[1][pos[0]][pos[1]][pos[2]] + curr_im[1][pos[0]][pos[1]][pos[2]-shift]);

            //for y
            m_Eng->volt[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(1,pos[0],pos[1],pos[2]);
            m_Eng->volt[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(1,pos[0],pos[1],pos[2])* ( m_Eng->curr[0][pos[0]][pos[1]][pos[2]] - m_Eng->curr[0][pos[0]][pos[1]][pos[2]-shift] - m_Eng->curr[2][pos[0]][pos[1]][pos[2]] + curr_outside[2]);
            volt_im[1][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(1,pos[0],pos[1],pos[2]);
            volt_im[1][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(1,pos[0],pos[1],pos[2])* ( curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im[0][pos[0]][pos[1]][pos[2]-shift] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im_outside[2]);

            //for z
            m_Eng->volt[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(2,pos[0],pos[1],pos[2]);
            m_Eng->volt[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(2,pos[0],pos[1],pos[2])* ( m_Eng->curr[1][pos[0]][pos[1]][pos[2]] - curr_outside[1] - m_Eng->curr[0][pos[0]][pos[1]][pos[2]] + m_Eng->curr[0][pos[0]][pos[1]-1][pos[2]]);
            volt_im[2][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(2,pos[0],pos[1],pos[2]);
            volt_im[2][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(2,pos[0],pos[1],pos[2])* ( curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im_outside[1]- curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im[0][pos[0]][pos[1]-1][pos[2]]);
    }
        pos[1] = 0;
        // moved the case y = pos[1] = 0 outside the loop
        // now treat the XZ-plane at the minimum x-coordinate and get the missing currents from pos[0] = maxX
        for (pos[0]=1; pos[0]<m_numLines[0]; ++pos[0]){
            curr_im_outside[0] = cos_ky*curr_im[0][pos[0]][maxY][pos[2]] - sin_ky*m_Eng->GetCurr(0,pos[0],maxY,pos[2]);
            curr_im_outside[2] = cos_ky*curr_im[2][pos[0]][maxY][pos[2]] - sin_ky*m_Eng->GetCurr(2,pos[0],maxY,pos[2]);
            curr_outside[0] =    cos_ky*curr_im[0][pos[0]][maxY][pos[2]] + sin_ky*m_Eng->GetCurr(0,pos[0],maxY,pos[2]);
            curr_outside[2] =    cos_ky*curr_im[2][pos[0]][maxY][pos[2]] + sin_ky*m_Eng->GetCurr(2,pos[0],maxY,pos[2]);

            m_Eng->volt[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(0,pos[0],pos[1],pos[2]);
            m_Eng->volt[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(0,pos[0],pos[1],pos[2])* ( m_Eng->curr[2][pos[0]][pos[1]][pos[2]] - curr_outside[2] - m_Eng->curr[1][pos[0]][pos[1]][pos[2]] + m_Eng->curr[1][pos[0]][pos[1]][pos[2]-shift]);
            volt_im[0][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(0,pos[0],pos[1],pos[2]);
            volt_im[0][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(0,pos[0],pos[1],pos[2])* ( curr_im[2][pos[0]][pos[1]][pos[2]] - curr_im_outside[2] - curr_im[1][pos[0]][pos[1]][pos[2]] + curr_im[1][pos[0]][pos[1]][pos[2]-shift]);
            //for y
            m_Eng->volt[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(1,pos[0],pos[1],pos[2]);
            m_Eng->volt[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(1,pos[0],pos[1],pos[2])* ( m_Eng->curr[0][pos[0]][pos[1]][pos[2]] - m_Eng->curr[0][pos[0]][pos[1]][pos[2]-shift] - m_Eng->curr[2][pos[0]][pos[1]][pos[2]] + m_Eng->curr[2][pos[0]-1][pos[1]][pos[2]]);
            volt_im[1][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(1,pos[0],pos[1],pos[2]);
            volt_im[1][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(1,pos[0],pos[1],pos[2])* ( curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im[0][pos[0]][pos[1]][pos[2]-shift] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im[2][pos[0]-1][pos[1]][pos[2]]);
            //for z
            m_Eng->volt[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(2,pos[0],pos[1],pos[2]);
            m_Eng->volt[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(2,pos[0],pos[1],pos[2])* ( m_Eng->curr[1][pos[0]][pos[1]][pos[2]] - m_Eng->curr[1][pos[0]-1][pos[1]][pos[2]] - m_Eng->curr[0][pos[0]][pos[1]][pos[2]] + curr_outside[0]);
            volt_im[2][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(2,pos[0],pos[1],pos[2]);
            volt_im[2][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(2,pos[0],pos[1],pos[2])* ( curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im[1][pos[0]-1][pos[1]][pos[2]]- curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im_outside[0]);
                    }
        pos[1] = 0;
        pos[0] = 0;
        // for x updates
        curr_outside[2] =    cos_ky*curr_im[2][pos[0]][maxY][pos[2]] + sin_ky*m_Eng->GetCurr(2,pos[0],maxY,pos[2]);
        curr_im_outside[2] = cos_ky*curr_im[2][pos[0]][maxY][pos[2]] - sin_ky*m_Eng->GetCurr(2,pos[0],maxY,pos[2]);
        m_Eng->volt[0][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(0,pos[0],pos[1],pos[2]);
        m_Eng->volt[0][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(0,pos[0],pos[1],pos[2])* ( m_Eng->curr[2][pos[0]][pos[1]][pos[2]] - curr_outside[2] - m_Eng->curr[1][pos[0]][pos[1]][pos[2]] + m_Eng->curr[1][pos[0]][pos[1]][pos[2]-shift]);
        volt_im[0][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(0,pos[0],pos[1],pos[2]);
        volt_im[0][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(0,pos[0],pos[1],pos[2])* ( curr_im[2][pos[0]][pos[1]][pos[2]] - curr_im_outside[2] - curr_im[1][pos[0]][pos[1]][pos[2]] + curr_im[1][pos[0]][pos[1]][pos[2]-shift]);
        // for y updates
        curr_outside[2] =    cos_kx*curr_im[2][maxX][pos[1]][pos[2]] + sin_kx*m_Eng->GetCurr(2,maxX,pos[1],pos[2]);
        curr_im_outside[2] = cos_kx*curr_im[2][maxX][pos[1]][pos[2]] - sin_kx*m_Eng->GetCurr(2,maxX,pos[1],pos[2]);
        m_Eng->volt[1][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(1,pos[0],pos[1],pos[2]);
        m_Eng->volt[1][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(1,pos[0],pos[1],pos[2])* ( m_Eng->curr[0][pos[0]][pos[1]][pos[2]] - m_Eng->curr[0][pos[0]][pos[1]][pos[2]-shift] - m_Eng->curr[2][pos[0]][pos[1]][pos[2]] + curr_outside[2]);
        volt_im[1][pos[0]][pos[1]][pos[2]] *=     m_Op_Pbc->GetVVedge(1,pos[0],pos[1],pos[2]);
        volt_im[1][pos[0]][pos[1]][pos[2]] +=     m_Op_Pbc->GetVIedge(1,pos[0],pos[1],pos[2])* ( curr_im[0][pos[0]][pos[1]][pos[2]] - curr_im[0][pos[0]][pos[1]][pos[2]-shift] - curr_im[2][pos[0]][pos[1]][pos[2]] + curr_im_outside[2]);

        // for z updates
        curr_outside[0] =    cos_ky*curr_im[0][pos[0]][maxY][pos[2]] + sin_ky*m_Eng->GetCurr(0,pos[0],maxY,pos[2]);
        curr_outside[1] =    cos_kx*curr_im[1][maxX][pos[1]][pos[2]] + sin_kx*m_Eng->GetCurr(1,maxX,pos[1],pos[2]);
        curr_im_outside[0] = cos_ky*curr_im[0][pos[0]][maxY][pos[2]] - sin_ky*m_Eng->GetCurr(0,pos[0],maxY,pos[2]);
        curr_im_outside[1] = cos_kx*curr_im[1][maxX][pos[1]][pos[2]] - sin_kx*m_Eng->GetCurr(1,maxX,pos[1],pos[2]);
        m_Eng->volt[2][pos[0]][pos[1]][pos[2]] *= m_Op_Pbc->GetVVedge(2,pos[0],pos[1],pos[2]);
        m_Eng->volt[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetVIedge(2,pos[0],pos[1],pos[2])* ( m_Eng->curr[1][pos[0]][pos[1]][pos[2]] - curr_outside[1] - m_Eng->curr[0][pos[0]][pos[1]][pos[2]] + curr_outside[0]);
        volt_im[2][pos[0]][pos[1]][pos[2]]     *= m_Op_Pbc->GetVVedge(2,pos[0],pos[1],pos[2]);
        volt_im[2][pos[0]][pos[1]][pos[2]]     += m_Op_Pbc->GetVIedge(2,pos[0],pos[1],pos[2])* ( curr_im[1][pos[0]][pos[1]][pos[2]] - curr_im_outside[1] - curr_im[0][pos[0]][pos[1]][pos[2]] + curr_im_outside[0]);
    }
};


void Engine_Ext_Pbc::DoPreVoltageUpdates(int threadID)
{
    DoVoltPhaseUpdates();

    // use periodic boundary condition to update the voltages at the lower boundary of the PBC directions
    // currents outside the FDTD domain are needed, but can be obtained from the currents on the opposite
    // end of the FDTD domain by multiplying them with a phase factor
    // i(nx=-1) = i(nx=Nx) * exp(-i k_x * Lx)
    // Re(i(-1)) = Re(i(Nx)) * cos(k_x * Lx) + Im(i(Nx)) * sin(k_x * Lx)
    // Im(i(-1)) = Im(i(Nx)) * cos(k_x * Lx) - Re(i(Nx)) * sin(k_x * Lx)

};
void Engine_Ext_Pbc::DoPreCurrentUpdates(int threadID)
{
    // use periodic boundary condition to update the currents at the upper boundary of the PBC directions
    // voltages outside the FDTD domain are needed, but can be obtained from the voltages on the opposite
    // end of the FDTD domain by multiplying them with a phase factor
    // u(nx=Nx) = u(nx=0) * exp(i k_x * Lx)
    // Re(u(Nx+1)) = Re(u(0)) * cos(k_x * Lx) - Im(u(u)) * sin(k_x * Lx)
    // Im(u(Nx+1)) = Im(u(0)) * cos(k_x * Lx) + Re(u(0)) * sin(k_x * Lx)
    //soft current excitation here (E-field excite)
    DoCurrPhaseUpdates();
};

void Engine_Ext_Pbc::Apply2Voltages(){
    //soft voltage excitation here (E-field excite)
    int exc_pos;
    unsigned int ny;
    unsigned int pos[3];
    int numTS = m_Eng->GetNumberOfTimesteps();
    unsigned int length = m_Op_Pbc->m_Exc->GetLength();
    FDTD_FLOAT* exc_volt_sinT =  m_Op_Pbc->m_Exc->GetVoltageSignal_s();
    FDTD_FLOAT* exc_volt_cosT =  m_Op_Pbc->m_Exc->GetVoltageSignal();

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
                m_Eng->Engine::SetVolt(ny,pos, m_Eng->Engine::GetVolt(ny,pos) + m_Op_Pbc->Volt_amp_cos[n]*exc_volt_cosT[exc_pos]/
                                                                              + m_Op_Pbc->Volt_amp_sin[n]*exc_volt_sinT[exc_pos]);
                volt_im[ny][pos[0]][pos[1]][pos[2]] += +m_Op_Pbc->Volt_amp_sin[n]*exc_volt_cosT[exc_pos] /
                                                       -m_Op_Pbc->Volt_amp_cos[n]*exc_volt_sinT[exc_pos];
            }
            break;
        }
     default:
    {
        break;
    }
    }
}

void Engine_Ext_Pbc::Apply2Current()
{
    int exc_pos;
    unsigned int ny;
    unsigned int pos[3];
    int numTS = m_Eng->GetNumberOfTimesteps();
    unsigned int length = m_Op_Pbc->m_Exc->GetLength();
    FDTD_FLOAT* exc_curr_sinT =  m_Op_Pbc->m_Exc->GetCurrentSignal_s();
    FDTD_FLOAT* exc_curr_cosT =  m_Op_Pbc->m_Exc->GetCurrentSignal(); // the signals time-dependence
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
                pos[2]=m_Op_Pbc->Curr_index[2][n]; // curr/volt_amp hold the spatial dependence, ex_curr/volt_... hold the temporal dependence
                m_Eng->Engine::SetCurr(ny,pos, m_Eng->Engine::GetCurr(ny,pos)  + m_Op_Pbc->Curr_amp_sin[n]*exc_curr_sinT[exc_pos]
                                                                               + m_Op_Pbc->Curr_amp_cos[n]*exc_curr_cosT[exc_pos]);
                curr_im[ny][pos[0]][pos[1]][pos[2]] += -m_Op_Pbc->Curr_amp_cos[n]*exc_curr_sinT[exc_pos] /
                                                       +m_Op_Pbc->Curr_amp_sin[n]*exc_curr_cosT[exc_pos];
            }
            break;
        }
     default:
    {
        break;
    }
    }
};


void Engine_Ext_Pbc::DoPostCurrentUpdates(int threadID){

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
                curr_im[2][pos[0]][pos[1]][pos[2]] += m_Op_Pbc->GetIV(2,pos[0],pos[1],pos[2]) * (volt_im[1][pos[0]][pos[1]][pos[2]] - volt_im[1][pos[0]+1][pos[1]][pos[2]] - volt_im[0][pos[0]][pos[1]][pos[2]] + volt_im[0][pos[0]][pos[1]+1][pos[2]]);
            }
        }
    }

};
void Engine_Ext_Pbc::DoPostVoltageUpdates(int threadID){
    for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
    {
        shift[0] = pos[0];
        for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
        {
            shift[1] = pos[1];
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
    }
};

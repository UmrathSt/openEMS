/* engine extension for periodic boundary conditions
*/
#include "engine_ext_pbc.h"
#include "operator_extension.h"

#include "FDTD/engine.h"

Engine_Ext_Pbc::Engine_Ext_Pbc(Operator_Ext_Pbc* op_ext) : Engine_Extension(op_ext)
{
    m_Op_Pbc = op_ext;
    m_numLines[0] = m_Op_Pbc->m_numLines[0];
    m_numLines[1] = m_Op_Pbc->m_numLines[1];
    m_ny = m_Op_Pbc->m_ny;
    m_nyP = m_Op_Pbc->m_nyP;
    m_nyPP = m_Op_Pbc->m_nyPP;
    m_LineNr = m_Op_Pbc->m_LineNr;
    m_LineNr_Shift = m_Op_Pbc->m_LineNr_Shift;

    m_Pbc_Coeff_nyP = m_Op_Pbc->m_Pbc_Coeff_nyP;
    m_Pbc_Coeff_nyPP = m_Op_Pbc->m_Pbc_Coeff_nyPP;

    m_volt_nyP = Create2DArray<FDTD_FLOAT>(m_numLines);
    m_volt_nyPP = Create2DArray<FDTD_FLOAT>(m_numLines);

    // ignore if some excitation is on the pbc
    // region
    SetNumberOfThreads(1);
}

Engine_Ext_Pbc::~Engine_Ext_Pbc()
{
    Delete2DArray(m_volt_nyP,m_numLines);
    m_volt_nyP = NULL;
    Delete2DArray(m_volt_nyPP,m_numLines);
    m_volt_nyPP = NULL;
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

void Engine_Ext_Pbc::DoPreVoltageUpdates(int threadID)
{
    if (IsActive()==false) return;
    if (m_Eng==NULL) return;
    if (threadID>=m_NrThreads)
        return;
    unsigned int pos[] = {0,0,0};
    unsigned int pos_shift[] = {0,0,0};
    pos[m_ny] = m_LineNr;
    pos_shift[m_ny] = m_LineNr_Shift;

    //switch for different engine types to access faster inline engine functions
    switch (m_Eng->GetType())
    {
    case Engine::BASIC:
        {
            for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
            {
                pos[m_nyP]=lineX+m_start.at(threadID);
                pos_shift[m_nyP] = pos[m_nyP];
                for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
                {
                    pos_shift[m_nyPP] = pos[m_nyPP];
                    m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] = m_Eng->Engine::GetVolt(m_nyP,pos_shift) - m_Op_Pbc->m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->Engine::GetVolt(m_nyP,pos);
                    m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] = m_Eng->Engine::GetVolt(m_nyPP,pos_shift) - m_Op_Pbc->m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->Engine::GetVolt(m_nyPP,pos);
                }
            }
            break;
        }
    case Engine::SSE:
        {
            Engine_sse* eng_sse = (Engine_sse*) m_Eng;
            for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
            {
                pos[m_nyP]=lineX+m_start.at(threadID);
                pos_shift[m_nyP] = pos[m_nyP];
                for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
                {
                    pos_shift[m_nyPP] = pos[m_nyPP];
                    m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] = eng_sse->Engine_sse::GetVolt(m_nyP,pos_shift) - m_Op_Pbc->m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * eng_sse->Engine_sse::GetVolt(m_nyP,pos);
                    m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] = eng_sse->Engine_sse::GetVolt(m_nyPP,pos_shift) - m_Op_Pbc->m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * eng_sse->Engine_sse::GetVolt(m_nyPP,pos);
                }
            }
            break;
        }
    default:
        for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
        {
            pos[m_nyP]=lineX+m_start.at(threadID);
            pos_shift[m_nyP] = pos[m_nyP];
            for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
            {
                pos_shift[m_nyPP] = pos[m_nyPP];
                m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] = m_Eng->GetVolt(m_nyP,pos_shift) - m_Op_Pbc->m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyP,pos);
                m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] = m_Eng->GetVolt(m_nyPP,pos_shift) - m_Op_Pbc->m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyPP,pos);
            }
        }
        break;
    }
}

void Engine_Ext_Pbc::DoPostVoltageUpdates(int threadID)
{
    if (IsActive()==false) return;
    if (m_Eng==NULL) return;
    if (threadID>=m_NrThreads)
        return;
    unsigned int pos[] = {0,0,0};
    unsigned int pos_shift[] = {0,0,0};
    pos[m_ny] = m_LineNr;
    pos_shift[m_ny] = m_LineNr_Shift;

    //switch for different engine types to access faster inline engine functions
    switch (m_Eng->GetType())
    {
    case Engine::BASIC:
        {
            for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
            {
                pos[m_nyP]=lineX+m_start.at(threadID);
                pos_shift[m_nyP] = pos[m_nyP];
                for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
                {
                    pos_shift[m_nyPP] = pos[m_nyPP];
                    m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] += m_Op_Pbc->m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->Engine::GetVolt(m_nyP,pos_shift);
                    m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] += m_Op_Pbc->m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->Engine::GetVolt(m_nyPP,pos_shift);
                }
            }
            break;
        }

    case Engine::SSE:
        {
            Engine_sse* eng_sse = (Engine_sse*) m_Eng;
            for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
            {
                pos[m_nyP]=lineX+m_start.at(threadID);
                pos_shift[m_nyP] = pos[m_nyP];
                for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
                {
                    pos_shift[m_nyPP] = pos[m_nyPP];
                    m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] += m_Op_Pbc->m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * eng_sse->Engine_sse::GetVolt(m_nyP,pos_shift);
                    m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] += m_Op_Pbc->m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * eng_sse->Engine_sse::GetVolt(m_nyPP,pos_shift);
                }
            }
            break;
        }

    default:
        for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
        {
            pos[m_nyP]=lineX+m_start.at(threadID);
            pos_shift[m_nyP] = pos[m_nyP];
            for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
            {
                pos_shift[m_nyPP] = pos[m_nyPP];
                m_volt_nyP[pos[m_nyP]][pos[m_nyPP]] += m_Op_Pbc->m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyP,pos_shift);
                m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]] += m_Op_Pbc->m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] * m_Eng->GetVolt(m_nyPP,pos_shift);
            }
        }
        break;
    }
}
void Engine_Ext_Pbc::Apply2Voltages(int threadID)
{
    if (IsActive()==false) return;
    if (threadID>=m_NrThreads)
        return;
    if (m_Eng==NULL) return;
    unsigned int pos[] = {0,0,0};
    pos[m_ny] = m_LineNr;

    //switch for different engine types to access faster inline engine functions
    switch (m_Eng->GetType())
    {
    case Engine::BASIC:
        {
            for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
            {
                pos[m_nyP]=lineX+m_start.at(threadID);
                for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
                {
                    m_Eng->Engine::SetVolt(m_nyP,pos, m_volt_nyP[pos[m_nyP]][pos[m_nyPP]]);
                    m_Eng->Engine::SetVolt(m_nyPP,pos, m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]]);
                }
            }
            break;
        }

    case Engine::SSE:
        {
            Engine_sse* eng_sse = (Engine_sse*) m_Eng;
            for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
            {
                pos[m_nyP]=lineX+m_start.at(threadID);
                for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
                {
                    eng_sse->Engine_sse::SetVolt(m_nyP,pos, m_volt_nyP[pos[m_nyP]][pos[m_nyPP]]);
                    eng_sse->Engine_sse::SetVolt(m_nyPP,pos, m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]]);
                }
            }
            break;
        }

    default:
        for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
        {
            pos[m_nyP]=lineX+m_start.at(threadID);
            for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
            {
                m_Eng->SetVolt(m_nyP,pos, m_volt_nyP[pos[m_nyP]][pos[m_nyPP]]);
                m_Eng->SetVolt(m_nyPP,pos, m_volt_nyPP[pos[m_nyP]][pos[m_nyPP]]);
            }
        }
        break;
    }

}

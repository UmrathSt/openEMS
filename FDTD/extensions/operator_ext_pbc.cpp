/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "operator_ext_pbc.h"
#include "engine_ext_pbc.h"

#include "tools/array_ops.h"
#include "CSPropMaterial.h"

Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op) : Operator_Extension(op)
{
   Initialize();
   for (unsigned int i  = 0; i<3; ++i){
       pbc_dirs[2*i] = false;
       pbc_dirs[2*i+1] = false;
       if (i==0 ){
           pbc_dirs[2*i] = true;
           pbc_dirs[2*i+1] = true;
       }
   }
   apply_PBC_to_operator(pbc_dirs);
}
Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext) : Operator_Extension(op, op_ext)
{
    Initialize();
    for (unsigned int i  = 0; i<3; ++i){
        pbc_dirs[2*i] = false;
        pbc_dirs[2*i+1] = false;
        if (i==0 )
            pbc_dirs[2*i] = true;
            pbc_dirs[2*i+1] = true;
    }
    apply_PBC_to_operator(pbc_dirs);
}
Operator_Ext_Pbc::~Operator_Ext_Pbc(){}

void Operator_Ext_Pbc::Initialize()
{
    m_numLines[0]= m_Op->GetNumberOfLines(0);
    m_numLines[1]= m_Op->GetNumberOfLines(1);
    m_numLines[2]= m_Op->GetNumberOfLines(2);
}
void Operator_Ext_Pbc::apply_PBC_to_operator(bool *dirs)
{
    for (int i=0; i<3; ++i){

    }

    cout << endl;
    unsigned int dir_lines[2] = {0, m_numLines[m_ny]-1};
    unsigned int pos[3];
    for(int i=0; i<3; ++i)
    {
        m_ny   = dirs[i];
        m_nyP  = (dirs[i]+1)%3;
        m_nyPP = (dirs[i]+2)%3;

        for (pos[m_nyP]=0; pos[m_nyP]<m_numLines[m_nyP]-1; ++pos[m_nyP])
        {

            for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[m_nyPP]-1; ++pos[m_nyPP])
            {
                if(dirs[2*i]){ // lowest mesh-line in direction i = (0,1,2) = ("x","y","z")
                    //cout << "operator_ext_pbc.cpp: setting the lower operator values according to PBC" << endl;
                    pos[m_ny] = 0;
                    m_Op->SetVV(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVV(m_nyP, pos[0], pos[1], pos[2],0);
                    m_Op->SetVV(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetII(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetII(m_nyP, pos[0], pos[1], pos[2],0);
                    m_Op->SetII(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyPP, pos[0], pos[1], pos[2], 0);

                }
                if(dirs[2*i+1]){ // lowest mesh-line in direction i = (0,1,2) = ("x","y","z")
                    //cout << "operator_ext_pbc.cpp: setting the higher operator values according to PBC" << endl;

                    pos[m_ny] = m_numLines[dirs[i]]-1;
                    m_Op->SetVV(m_ny, pos[0], pos[1], pos[2], 1);
                    m_Op->SetVV(m_nyP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetVV(m_nyPP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetVI(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetII(m_ny, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_nyP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_nyPP, pos[0], pos[1], pos[2], 1);
                    pos[m_ny] -= 1; // outside the FDTD domain...?
                    m_Op->SetIV(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetII(m_ny, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_nyP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_nyPP, pos[0], pos[1], pos[2], 1);


                }

            }
        }
    }
}




Engine_Extension* Operator_Ext_Pbc::CreateEngineExtention()
{
    cout << "operator_ext_pbc.cpp: created engine extenesion engine_ext_pbc(this)" << endl;
    Engine_Ext_Pbc* eng_ext_pbc = new Engine_Ext_Pbc(this);
    return eng_ext_pbc;
}
Operator_Extension* Operator_Ext_Pbc::Clone(Operator* op)
{
    if (dynamic_cast<Operator_Ext_Pbc*>(this)==NULL)
        return NULL;
    return new Operator_Ext_Pbc(op, this);
}
void Operator_Ext_Pbc::Set_k_pbc(FDTD_FLOAT *k){
    for (int i=0; i<3; ++i){
        k_pbc[i] = k[i];
        cout << "operator_ext_pbc.cpp: now k_pbc[i] is equal to = " << k_pbc[i] << endl;
    }
}

bool Operator_Ext_Pbc::BuildExtension()
{
    m_Exc = m_Op->GetExcitationSignal();
    unsigned int m_numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

    if (m_Op->k_pbc[0] == 0 && m_Op->k_pbc[1] == 0 && m_Op->k_pbc[2] == 0)
    {
        cerr << "Operator_Ext_Pbc::BuildExtension: Warning, Obviously the PBC-Extension was used without setting m_k_PBC, trying to proceed anyways ..." << endl;
        return false;
    }
    return true;
}

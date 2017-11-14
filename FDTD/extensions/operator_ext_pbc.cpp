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
       if (i==0){
           apply_PBC_to_operator(i);
       }
   }
}
Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext) : Operator_Extension(op, op_ext)
{
    Initialize();
    for (unsigned int i  = 0; i<3; ++i){
        if (i==0 || i==1)
            apply_PBC_to_operator(i);
    }
}
Operator_Ext_Pbc::~Operator_Ext_Pbc(){}

void Operator_Ext_Pbc::Initialize()
{
    m_numLines[0]= m_Op->GetNumberOfLines(0);
    m_numLines[1]= m_Op->GetNumberOfLines(1);
    m_numLines[2]= m_Op->GetNumberOfLines(2);
}
void Operator_Ext_Pbc::apply_PBC_to_operator(unsigned int dirs)
{
    for (int n=0; n<3; ++n){
        m_ny   = dirs;
        m_nyP  = (dirs+1)%3;
        m_nyPP = (dirs+2)%3;
        unsigned int dir_lines[2] = {0, m_numLines[m_ny]};
        unsigned int pos[3];
        // set operator values in the +- nyP-nyPP-plane to 1
        for (int i=0; i<2; ++i){
            for (pos[m_nyP]=0; pos[m_nyP]<m_numLines[m_nyP]; ++pos[m_nyP]){
                for (pos[m_nyPP]=0; pos[m_nyPP]>m_numLines[m_nyPP]; ++pos[m_nyPP]){
                    pos[m_ny] = dir_lines[i];
                    m_Op->SetVV(m_ny, pos[0], pos[1], pos[2], 1);
                    m_Op->SetVV(m_nyP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetVV(m_nyPP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_ny, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_nyP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetII(m_nyPP, pos[0], pos[1], pos[2], 1);
                    m_Op->SetVI(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetVI(m_nyPP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_ny, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyP, pos[0], pos[1], pos[2], 0);
                    m_Op->SetIV(m_nyPP, pos[0], pos[1], pos[2], 0);
                }
            }
        }
    }
}


Engine_Extension* Operator_Ext_Pbc::CreateEngineExtention()
{
    Engine_Ext_Pbc* eng_ext_pbc = new Engine_Ext_Pbc(this);
    return eng_ext_pbc;
}
Operator_Extension* Operator_Ext_Pbc::Clone(Operator* op)
{
    if (dynamic_cast<Operator_Ext_Pbc*>(this)==NULL)
        return NULL;
    return new Operator_Ext_Pbc(op, this);
}

bool Operator_Ext_Pbc::BuildExtension()
{
    unsigned int m_numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

    if (m_Op->k_PBC[0] == 0 && m_Op->k_PBC[1] == 0 && m_Op->k_PBC[2] == 0)
    {
        cerr << "Operator_Ext_Pbc::BuildExtension: Warning, Obviously the PBC-Extension was used without setting m_k_PBC, aborting ..." << endl;
        return false;
    }
    return true;
}

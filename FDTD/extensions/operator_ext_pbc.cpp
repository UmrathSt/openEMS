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
       if (m_Op->dir_is_pbc[i]){
           cout << "direction " << i << " is PBC" << endl;
           cout << "namely: k[" << i << "]= " << m_Op->k_PBC[i] << endl;
           apply_PBC_to_operator(i);
       }
   }
}
Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext) : Operator_Extension(op, op_ext)
{
    Initialize();
    for (unsigned int i  = 0; i<3; ++i){
        if (m_Op->dir_is_pbc[i])
            apply_PBC_to_operator(i);
    }
}
Operator_Ext_Pbc::~Operator_Ext_Pbc(){}

void Operator_Ext_Pbc::Initialize()
{
    m_numLines[0]= m_Op->GetNumberOfLines(0);
    m_numLines[1]= m_Op->GetNumberOfLines(1);
    m_numLines[2]= m_Op->GetNumberOfLines(2);
    for(int i = 0; i<3; ++i){
        if (m_Op->dir_is_pbc[i]){
            cout << "k_PBC[" << i << "] = " << m_Op->k_PBC[i] << endl;
        }
    }
}
void Operator_Ext_Pbc::apply_PBC_to_operator(unsigned int dir)
{
    m_ny   = dir;
    m_nyP  = (dir+1)%3;
    m_nyPP = (dir+2)%3;
    unsigned int dir_lines[2] = {0, m_numLines[m_ny]-1};
    unsigned int pos[3];
    // set operator values in the +- nyP-nyPP-plane to 1
    for (int i=0; i<2; ++i ){
        pos[m_ny] = dir_lines[i];
        for (pos[m_nyP]=0; pos[m_nyP]<m_numLines[m_nyP]-1; ++pos[m_nyP]){
            for (pos[m_nyPP]=0; pos[m_nyPP]>m_numLines[m_nyPP]-1; ++pos[m_nyPP]){
                m_Op->SetVV(0, pos[0], pos[1], pos[2], 1);
                m_Op->SetVV(1, pos[0], pos[1], pos[2], 1);
                m_Op->SetVV(2, pos[0], pos[1], pos[2], 1);
                m_Op->SetII(0, pos[0], pos[1], pos[2], 1);
                m_Op->SetII(1, pos[0], pos[1], pos[2], 1);
                m_Op->SetII(2, pos[0], pos[1], pos[2], 1);
            }

        }

    }
    cout << "PBC's successfully applied to Operator!" << endl;
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
        exit(-3);
        return false;
    }
    return true;
}

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
}
Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext) : Operator_Extension(op, op_ext)
{
    Initialize();
    setKParallel(op_ext->kparallel);

}
Operator_Ext_Pbc::~Operator_Ext_Pbc(){}

void Operator_Ext_Pbc::Initialize()
{

    m_numLines[0]= m_Op->GetNumberOfLines(0);
    m_numLines[1]= m_Op->GetNumberOfLines(1);
    m_numLines[2]= m_Op->GetNumberOfLines(2);
    for(int i = 0; i<3; ++i){
        kparallel[i] = -1;
    }
}

void Operator_Ext_Pbc::setKParallel(float *kpar)
{
    kparallel[0] = kpar[0];
    kparallel[1] = kpar[1];
    kparallel[2] = kpar[2];
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
    for(int i=0; i<3; ++i)
        kparallel[i] = m_Op->m_k_PBC[i];
    unsigned int m_numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

    if (kparallel[0] == -1 && kparallel[1] == -1 && kparallel[2] == -1)
    {
        cerr << "Operator_Ext_Pbc::BuildExtension: Warning, Extension not initialized! Use SetKParallel!! Abort build!!" << endl;
        return false;
    }
    return true;
}

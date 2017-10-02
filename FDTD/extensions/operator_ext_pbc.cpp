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

Operator_Ext_Pbc::~Operator_Ext_Pbc()
{
}

//bool Operator_Ext_Pbc::IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const
//{
//    if ((m_ny==0) && (!m_top) && (R0_included || closedAlpha))
//        return false;
//    if ((m_ny==1) && (closedAlpha))
//        return false;
//    return true; // to be modiefied, does not yet make sense
//}

//bool Operator_Ext_Pbc::IsCylindricalMultiGridSave(bool child) const
//{
//    if (m_ny==2) //always allow in z-direction
//        return true;
//    if ((m_ny==0) && (m_top) && (!child)) //if top r-direction and is not a child grid allow Mur...
//        return true;
//    //in all other cases this ABC is not save to use in CylindricalMultiGrid
//    return false; // to be modified, does not yet make sense
//}

Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op) : Operator_Extension(op)
{
    Initialize();
}

void Operator_Ext_Pbc::Initialize()
{
    kparallel = {-1, -1, -1};
}

void Operator_Ext_Pbc::SetKParallel(std::vector<double> &kpar)
{
    kparallel = kpar;
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
    if (kparallel[0] == -1 && kparallel[1] == -1 && kparallel[2] == -1)
    {
        cerr << "Operator_Ext_Pbc::BuildExtension: Warning, Extension not initialized! Use SetKParallel!! Abort build!!" << endl;
        return false;
    }
    double dT = m_Op->GetTimestep();
    unsigned int pos[] = {0,0,0};
    pos[m_ny] = m_LineNr;
    double delta = fabs(m_Op->GetEdgeLength(m_ny,pos));
    double coord[] = {0,0,0};
    coord[0] = m_Op->GetDiscLine(0,pos[0]);
    coord[1] = m_Op->GetDiscLine(1,pos[1]);
    coord[2] = m_Op->GetDiscLine(2,pos[2]);

    double eps,mue;
    double c0t;

    if (m_LineNr==0)
        coord[m_ny] = m_Op->GetDiscLine(m_ny,pos[m_ny]) + delta/2 / m_Op->GetGridDelta();
    else
        coord[m_ny] = m_Op->GetDiscLine(m_ny,pos[m_ny]) - delta/2 / m_Op->GetGridDelta();

    int posBB[3];
    posBB[m_ny]  =pos[m_ny];
    posBB[m_nyPP]=-1;

    for (pos[m_nyP]=0; pos[m_nyP]<m_numLines[0]; ++pos[m_nyP])
    {
        posBB[m_nyP]=pos[m_nyP];
        vector<CSPrimitives*> vPrims = m_Op->GetPrimitivesBoundBox(posBB[0], posBB[1], posBB[2], CSProperties::MATERIAL);
        coord[m_nyP] = m_Op->GetDiscLine(m_nyP,pos[m_nyP]);
        for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[1]; ++pos[m_nyPP])
        {
            coord[m_nyPP] = m_Op->GetDiscLine(m_nyPP,pos[m_nyPP]);
//			CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, CSProperties::MATERIAL, false);
            CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, false);
            if (prop)
            {
                CSPropMaterial* mat = prop->ToMaterial();

                //nP
                eps = mat->GetEpsilonWeighted(m_nyP,coord);
                mue = mat->GetMueWeighted(m_nyP,coord);
                if (m_v_phase>0.0)
                    c0t = m_v_phase * dT;
                else
                    c0t = __C0__ * dT / sqrt(eps*mue);
                m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] = (c0t - delta) / (c0t + delta);

                //nPP
                eps = mat->GetEpsilonWeighted(m_nyPP,coord);
                mue = mat->GetMueWeighted(m_nyPP,coord);
                if (m_v_phase>0.0)
                    c0t = m_v_phase * dT;
                else
                    c0t = __C0__ * dT / sqrt(eps*mue);
                m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] = (c0t - delta) / (c0t + delta);

            }
            else
            {
                if (m_v_phase>0.0)
                    c0t = m_v_phase * dT;
                else
                    c0t = __C0__ / sqrt(m_Op->GetBackgroundEpsR()*m_Op->GetBackgroundMueR()) * dT;
                m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] = (c0t - delta) / (c0t + delta);
                m_Pbc_Coeff_nyPP[pos[m_nyP]][pos[m_nyPP]] = m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]];
            }
//			cerr << m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] << " : " << m_Pbc_Coeff_nyP[pos[m_nyP]][pos[m_nyPP]] << endl;
        }
    }
//	cerr << "Operator_Ext_Pbc_ABC::BuildExtension(): " << m_ny << " @ " << m_LineNr << endl;
    return true;
}

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
#include "ContinuousStructure.h"
#include "CSPrimCurve.h"
#include "CSPropPBCExcitation.h"
#include "CSPropExcitation.h"

Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op) : Operator_Extension(op)
{
    Init();
    copy_operator_vals();
    apply_PBC_to_operator(pbc_dirs);
}
Operator_Ext_Pbc::Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext) : Operator_Extension(op, op_ext)
{
    Init();
    copy_operator_vals();
    apply_PBC_to_operator(pbc_dirs);
}
Operator_Ext_Pbc::~Operator_Ext_Pbc(){}

void Operator_Ext_Pbc::copy_operator_vals()
{
    VV = m_Op->vv;
    VI = m_Op->vi;
    IV = m_Op->iv;
    II = m_Op->ii;
};


void Operator_Ext_Pbc::Init()
{
    m_Exc = m_Op->GetExcitationSignal();
    m_numLines[0]= m_Op->GetNumberOfLines(0);
    m_numLines[1]= m_Op->GetNumberOfLines(1);
    m_numLines[2]= m_Op->GetNumberOfLines(2);
    VV = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
    VI = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
    IV = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
    II = Create_N_3DArray<FDTD_FLOAT>(m_numLines);
    std::cout << "I initialized the Operator_Ext_Pbc and created VV of (" << 3 << "," << m_numLines[0] << "," << m_numLines[1] << "," << m_numLines[2] << ")" << std::endl;
    Operator_Extension::Init();
    Volt_delay = 0;
    Volt_amp_sin = 0; // time dependence
    Volt_amp_cos = 0;
    Volt_dir = 0;
    Volt_Count = 0;
    Curr_delay = 0;
    Curr_amp_sin = 0;
    Curr_amp_cos = 0;
    Curr_dir = 0;
    Curr_Count = 0;

    for (int n=0; n<3; ++n)
    {
        Volt_index[n] = 0;
        Curr_index[n] = 0;
        Volt_Count_Dir[n] = 0;
        Curr_Count_Dir[n] = 0;
    }
    BuildExtension();
}

void Operator_Ext_Pbc::Reset()
{
    Operator_Extension::Reset();
    delete[] Volt_delay;
    Volt_delay = 0;
    delete[] Volt_dir;
    Volt_dir = 0;
    delete[] Volt_amp_sin;
    Volt_amp_sin = 0;
    delete[] Volt_amp_cos;
    Volt_amp_cos = 0;
    delete[] Curr_delay;
    Curr_delay = 0;
    delete[] Curr_dir;
    Curr_dir = 0;
    delete[] Curr_amp_sin;
    Curr_amp_sin = 0;
    delete[] Curr_amp_cos;
    Curr_amp_cos = 0;

    Volt_Count = 0;
    Curr_Count = 0;

    for (int n=0; n<3; ++n)
    {
        delete[] Volt_index[n];
        Volt_index[n] = 0;
        delete[] Curr_index[n];
        Curr_index[n] = 0;

        Volt_Count_Dir[n] = 0;
        Curr_Count_Dir[n] = 0;
    }
}
void Operator_Ext_Pbc::apply_PBC_to_operator(bool *dirs)
{


    for(int i=0; i<3; ++i)
    {
        m_ny   = dirs[i];
        m_nyP  = (dirs[i]+1)%3;
        m_nyPP = (dirs[i]+2)%3;
        FDTD_FLOAT pp_val = 1;
        FDTD_FLOAT pq_val = 0;
        // set the operator to zero/one such that the updates at the boundary
        // can be treated separately in engine_ext_pbc.cpp
        for (pos[m_nyP]=0; pos[m_nyP]<m_numLines[m_nyP]; ++pos[m_nyP])
        {
            for (pos[m_nyPP]=0; pos[m_nyPP]<m_numLines[m_nyPP]; ++pos[m_nyPP])
            {
                if(dirs[2*i]){ // The voltages have missing current neighbours at the lower boundary
                    pos[m_ny] = 0;

                    m_Op->SetVV(m_ny, pos[0], pos[1], pos[2],   pp_val);
                    m_Op->SetVV(m_nyP, pos[0], pos[1], pos[2],  pp_val);
                    m_Op->SetVV(m_nyPP, pos[0], pos[1], pos[2], pp_val);
                    m_Op->SetVI(m_ny, pos[0], pos[1], pos[2],   pq_val);
                    m_Op->SetVI(m_nyP, pos[0], pos[1], pos[2],  pq_val);
                    m_Op->SetVI(m_nyPP, pos[0], pos[1], pos[2], pq_val);
                }
                 if(dirs[2*i+1]){ // The currents have missing voltage neighbours at the higher boundaries
                    pos[m_ny] = m_numLines[dirs[i]]-1; // and are therefore updated separately

                    m_Op->SetII(m_ny, pos[0], pos[1], pos[2],   pp_val);
                    m_Op->SetII(m_nyP, pos[0], pos[1], pos[2],  pp_val);
                    m_Op->SetII(m_nyPP, pos[0], pos[1], pos[2], pp_val);
                    m_Op->SetIV(m_ny, pos[0], pos[1], pos[2],   pq_val);
                    m_Op->SetIV(m_nyP, pos[0], pos[1], pos[2],  pq_val);
                    m_Op->SetIV(m_nyPP, pos[0], pos[1], pos[2], pq_val);
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
void Operator_Ext_Pbc::Set_k_pbc(FDTD_FLOAT *k){
    for (int i=0; i<3; ++i){
        k_pbc[i] = k[i];
    }
}

void Operator_Ext_Pbc::Set_pbc_dirs(bool *dirs){
    for (int i=0; i<6; ++i){
        pbc_dirs[i] = dirs[i];
    }
}

bool Operator_Ext_Pbc::BuildExtension()
{
    Build_PBCExcitation();

    unsigned int m_numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

    if (m_Op->k_pbc[0] == 0 && m_Op->k_pbc[1] == 0 && m_Op->k_pbc[2] == 0)
    {
        cerr << "Operator_Ext_Pbc::BuildExtension: Warning, Obviously the PBC-Extension was used without setting m_k_PBC, trying to proceed anyways ..." << endl;
        return false;
    }
    return true;
}

bool Operator_Ext_Pbc::Build_PBCExcitation()
{

    double dT = m_Op->GetTimestep();
    if (dT==0)
        return false;
    if (m_Exc==0)
        return false;

    Reset();
    ContinuousStructure* CSX = m_Op->GetGeometryCSX();

    unsigned int pos[3];
    double amp_cos=0;
    double amp_sin=0;
    vector<unsigned int> volt_vIndex[3];
    vector<FDTD_FLOAT> volt_vExcit_sin;
    vector<FDTD_FLOAT> volt_vExcit_cos;
    vector<unsigned int> volt_vDelay;
    vector<unsigned int> volt_vDir;
    double volt_coord[3];

    vector<unsigned int> curr_vIndex[3];
    vector<FDTD_FLOAT> curr_vExcit_sin;
    vector<FDTD_FLOAT> curr_vExcit_cos;
    vector<unsigned int> curr_vDelay;
    vector<unsigned int> curr_vDir;
    double curr_coord[3];

    vector<CSProperties*> vec_prop = CSX->GetPropertyByType(CSProperties::PBCEXCITATION);

    if (vec_prop.size()==0)
    {
        cerr << "Operator::CalcFieldExcitation: Warning, no PBC excitation properties found" << endl;
        return false;
    }

    CSPropPBCExcitation* elec=NULL;
    CSProperties* prop=NULL;
    int priority=0;
    unsigned int numLines[] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

    for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
    {
        for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
        {
            for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
            {
                //electric field excite
                for (int n=0; n<3; ++n)
                {

                    if (m_Op->GetYeeCoords(n,pos,volt_coord,false)==false)
                        continue;
                    if (m_CC_R0_included && (n==2) && (pos[0]==0))
                        volt_coord[1] = m_Op->GetDiscLine(1,0);

                    if (m_CC_R0_included && (n==1) && (pos[0]==0))
                        continue;
                    for (size_t p=0; p<vec_prop.size(); ++p)
                    {
                        prop = vec_prop.at(p);
                        elec = prop->ToPBCExcitation();
                        if (elec==NULL){
                            continue;
                        }
                        if (prop->CheckCoordInPrimitive(volt_coord,priority,true))
                        {                            
                            if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==0) || (elec->GetExcitType()==1) ))//&& (pos[n]<numLines[n]-1))
                            {
                                amp_sin = elec->GetWeightedExcitation(n,volt_coord,1)*m_Op->GetEdgeLength(n,pos);// delta[n]*gridDelta;
                                amp_cos = elec->GetWeightedExcitation(n,volt_coord,0)*m_Op->GetEdgeLength(n,pos);// delta[n]*gridDelta;
                                if (amp_sin!=0 || amp_cos!=0)
                                {
                                    volt_vExcit_sin.push_back(amp_sin);
                                    volt_vExcit_cos.push_back(amp_cos);
                                    volt_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
                                    volt_vDir.push_back(n);
                                    volt_vIndex[0].push_back(pos[0]);
                                    volt_vIndex[1].push_back(pos[1]);
                                    volt_vIndex[2].push_back(pos[2]);
                                }

                                if (elec->GetExcitType()==1) //hard excite
                                {
                                    m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
                                    m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
                                }
                            }
                    }
                    }
                }

                //magnetic field excite
                for (int n=0; n<3; ++n)
                {
                    if ((pos[0]>=numLines[0]-1) || (pos[1]>=numLines[1]-1) || (pos[2]>=numLines[2]-1))
                        continue;  //skip the last H-Line which is outside the FDTD-domain
                    if (m_Op->GetYeeCoords(n,pos,curr_coord,true)==false)
                        continue;
                    for (size_t p=0; p<vec_prop.size(); ++p)
                    {
                        prop = vec_prop.at(p);
                        elec = prop->ToPBCExcitation();
                        if (elec==NULL)
                            continue;
                        if (prop->CheckCoordInPrimitive(curr_coord,priority,true))
                        {
                            if ((elec->GetActiveDir(n)) && ( (elec->GetExcitType()==2) || (elec->GetExcitType()==3) ))
                            {
                                amp_sin = elec->GetWeightedExcitation(n,curr_coord,1)*m_Op->GetEdgeLength(n,pos,true);// delta[n]*gridDelta;
                                amp_cos = elec->GetWeightedExcitation(n,curr_coord,0)*m_Op->GetEdgeLength(n,pos,true);// delta[n]*gridDelta;



                                if (amp_cos!=0 || amp_sin!=0)
                                {
                                    curr_vExcit_sin.push_back(amp_sin);
                                    curr_vExcit_cos.push_back(amp_cos);
                                    curr_vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
                                    curr_vDir.push_back(n);
                                    curr_vIndex[0].push_back(pos[0]);
                                    curr_vIndex[1].push_back(pos[1]);
                                    curr_vIndex[2].push_back(pos[2]);
                                }
                                if (elec->GetExcitType()==3) //hard excite
                                {
                                    m_Op->SetII(n,pos[0],pos[1],pos[2], 0 );
                                    m_Op->SetIV(n,pos[0],pos[1],pos[2], 0 );
                                }
                            }
                        }
                    }
                }

            }
        }
    }




    // set voltage excitations
    setupVoltageExcitation( volt_vIndex, volt_vExcit_sin, volt_vExcit_cos, volt_vDelay, volt_vDir );

    // set current excitations
    setupCurrentExcitation( curr_vIndex, curr_vExcit_sin, curr_vExcit_cos, curr_vDelay, curr_vDir );

    return true;
}


void Operator_Ext_Pbc::setupVoltageExcitation( vector<unsigned int> const volt_vIndex[3], vector<FDTD_FLOAT> const& volt_vExcit_sin,
        vector<FDTD_FLOAT> const& volt_vExcit_cos, vector<unsigned int> const& volt_vDelay, vector<unsigned int> const& volt_vDir )
{
    Volt_Count = volt_vIndex[0].size();
    cout << "---------------------" << endl;
    cout << "Operator_ext_Pbc::setupVoltageExcitation: Volt_count = " << Volt_Count << endl;
    cout << "---------------------" << endl;
    for (int n=0; n<3; n++)
    {
        Volt_Count_Dir[n]=0;
        delete[] Volt_index[n];
        Volt_index[n] = new unsigned int[Volt_Count];
    }
    delete[] Volt_delay;
    delete[] Volt_amp_sin;
    delete[] Volt_amp_cos;
    delete[] Volt_dir;
    Volt_delay = new unsigned int[Volt_Count];
    Volt_amp_sin = new FDTD_FLOAT[Volt_Count];
    Volt_amp_cos = new FDTD_FLOAT[Volt_Count];
    Volt_dir = new unsigned short[Volt_Count];

//	cerr << "Excitation::setupVoltageExcitation(): Number of voltage excitation points: " << Volt_Count << endl;
//	if (Volt_Count==0)
//		cerr << "No E-Field/voltage excitation found!" << endl;
    for (int n=0; n<3; n++)
        for (unsigned int i=0; i<Volt_Count; i++)
            Volt_index[n][i] = volt_vIndex[n].at(i);
    for (unsigned int i=0; i<Volt_Count; i++)
    {
        Volt_delay[i] = volt_vDelay.at(i);
        Volt_amp_sin[i]   = volt_vExcit_sin.at(i);
        Volt_amp_cos[i]   = volt_vExcit_cos.at(i);
        Volt_dir[i]   = volt_vDir.at(i);
        ++Volt_Count_Dir[Volt_dir[i]];
    }
}

void Operator_Ext_Pbc::setupCurrentExcitation( vector<unsigned int> const curr_vIndex[3], vector<FDTD_FLOAT> const& curr_vExcit_sin,
        vector<FDTD_FLOAT> const& curr_vExcit_cos, vector<unsigned int> const& curr_vDelay, vector<unsigned int> const& curr_vDir )
{
    Curr_Count = curr_vIndex[0].size();
    for (int n=0; n<3; n++)
    {
        Curr_Count_Dir[n]=0;
        delete[] Curr_index[n];
        Curr_index[n] = new unsigned int[Curr_Count];
    }
    delete[] Curr_delay;
    delete[] Curr_amp_sin;
    delete[] Curr_amp_cos;
    delete[] Curr_dir;
    Curr_delay = new unsigned int[Curr_Count];
    Curr_amp_sin = new FDTD_FLOAT[Curr_Count];
    Curr_amp_cos = new FDTD_FLOAT[Curr_Count];
    Curr_dir = new unsigned short[Curr_Count];

//	cerr << "Excitation::setupCurrentExcitation(): Number of current excitation points: " << Curr_Count << endl;
//	if (Curr_Count==0)
//		cerr << "No H-Field/current excitation found!" << endl;
    for (int n=0; n<3; ++n)
        for (unsigned int i=0; i<Curr_Count; i++)
            Curr_index[n][i] = curr_vIndex[n].at(i);
    for (unsigned int i=0; i<Curr_Count; i++)
    {
        Curr_delay[i] = curr_vDelay.at(i);
        Curr_amp_sin[i]   = curr_vExcit_sin.at(i);
        Curr_amp_cos[i]   = curr_vExcit_cos.at(i);
        Curr_dir[i]   = curr_vDir.at(i);
        ++Curr_Count_Dir[Curr_dir[i]];
    }

}
void Operator_Ext_Pbc::ShowStat(ostream &ostr)  const
{
    Operator_Extension::ShowStat(ostr);
    cout << "Sin Voltage excitations\t: " << Volt_Count    << "\t (" << Volt_Count_Dir[0] << ", " << Volt_Count_Dir[1] << ", " << Volt_Count_Dir[2] << ")" << endl;
    cout << "Cos Voltage excitations\t: " << Volt_Count    << "\t (" << Volt_Count_Dir[0] << ", " << Volt_Count_Dir[1] << ", " << Volt_Count_Dir[2] << ")" << endl;
    cout << "Sin Current excitations\t: " << Curr_Count << "\t (" << Curr_Count_Dir[0] << ", " << Curr_Count_Dir[1] << ", " << Curr_Count_Dir[2] << ")" << endl;
    cout << "Cos Current excitations\t: " << Curr_Count << "\t (" << Curr_Count_Dir[0] << ", " << Curr_Count_Dir[1] << ", " << Curr_Count_Dir[2] << ")" << endl;
    cout << "Excitation Length (TS)\t: " << m_Exc->GetLength() << endl;
    cout << "Excitation Length (s)\t: " << m_Exc->GetLength()*m_Op->GetTimestep() << endl;
}

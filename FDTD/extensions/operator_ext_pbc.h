#ifndef OPERATOR_EXT_PBC_H
#define OPERATOR_EXT_PBC_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Operator_Ext_Pbc : public Operator_Extension
{
    friend class Engine_Ext_Pbc;
    friend class openEMS;
public:
    Operator_Ext_Pbc(Operator* op);
    ~Operator_Ext_Pbc();
    inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->GetVV(n,x,y,z); };
    inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->GetVI(n,x,y,z); };
    inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->GetIV(n,x,y,z); };
    inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->GetII(n,x,y,z); };
    inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT val) const { return m_Op->SetVV(n,x,y,z,val); };
    inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT val) const { return m_Op->SetVI(n,x,y,z,val); };
    inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT val) const { return m_Op->SetIV(n,x,y,z,val); };
    inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT val) const { return m_Op->SetII(n,x,y,z,val); };
    virtual Operator_Extension* Clone(Operator* op);
    // sets the phase difference between opposite sides of the periodic structure
    // Example: periodicity (exp(i * k * r) in x-direction only -> m_k_PBC = {1, 0, 0};
    void setupCurrentExcitation( vector<unsigned int> const curr_vIndex[3], vector<FDTD_FLOAT> const& curr_vExcit_sin,
            vector<FDTD_FLOAT> const& curr_vExcit_cos, vector<unsigned int> const& curr_vDelay, vector<unsigned int> const& curr_vDir );
    void setupVoltageExcitation( vector<unsigned int> const volt_vIndex[3], vector<FDTD_FLOAT> const& volt_vExcit,
            vector<FDTD_FLOAT> const& volt_vExcit_cos, vector<unsigned int> const& volt_vDelay, vector<unsigned int> const& volt_vDir );
    bool Build_PBCExcitation();
    void apply_PBC_to_operator(bool *dirs);
    void Set_k_pbc(FDTD_FLOAT *k_pbc);
    void Set_pbc_dirs(bool *dirs);

    virtual bool BuildExtension();
    virtual Engine_Extension* CreateEngineExtention();

    virtual bool IsMPISave() const {return true;}

    virtual string GetExtensionName() const {return string("PeriodicBoundaryCondition Extension");}

    virtual void ShowStat(ostream &ostr) const;

protected:
    Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext);
    void Init();
    void Reset();
    unsigned int m_numLines[3];
    unsigned int pos[3];
    int m_ny;
    int m_nyP,m_nyPP;
    bool pbc_dirs[6] = {0};
    FDTD_FLOAT k_pbc[3] = {0};
    bool pbc_used = true;

    Excitation* m_Exc;
    //E-Field/voltage Excitation
    unsigned int Volt_Count;
    unsigned int Volt_Count_Dir[3];
    unsigned int* Volt_index[3];
    unsigned short* Volt_dir;
    FDTD_FLOAT* Volt_amp_sin; //represented as edge-voltages!!
    FDTD_FLOAT* Volt_amp_cos; //represented as edge-voltages!!
    unsigned int* Volt_delay;

    //H-Field/current Excitation
    unsigned int Curr_Count;
    unsigned int Curr_Count_Dir[3];
    unsigned int* Curr_index[3];
    unsigned short* Curr_dir;
    FDTD_FLOAT* Curr_amp_cos; //represented as edge-currents!!
    FDTD_FLOAT* Curr_amp_sin; //represented as edge-currents!!
    unsigned int* Curr_delay;

};


#endif // OPERATOR_EXT_PBC_H

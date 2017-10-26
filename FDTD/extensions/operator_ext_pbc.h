#ifndef OPERATOR_EXT_PBC_H
#define OPERATOR_EXT_PBC_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Operator_Ext_Pbc : public Operator_Extension
{
    friend class Engine_Ext_Pbc;
public:
    Operator_Ext_Pbc(Operator* op);
    ~Operator_Ext_Pbc();
    inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->vv[n][x][y][z]; }
    inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->vi[n][x][y][z]; }
    inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->iv[n][x][y][z]; }
    inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { return m_Op->ii[n][x][y][z]; }
    virtual Operator_Extension* Clone(Operator* op);
    // sets the phase difference between opposite sides of the periodic structure
    // Example: periodicity (exp(i * k * r) in x-direction only -> kparallel = {1, 0, 0};
    void setKParallel(float *kpar);

    virtual bool BuildExtension();
    virtual Engine_Extension* CreateEngineExtention();

    virtual bool IsMPISave() const {return true;}

    virtual string GetExtensionName() const {return string("PeriodicBoundaryCondition Extension");}



protected:
    Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext);
    void Initialize();
    float kparallel[3];
    unsigned int m_numLines[3];
};


#endif // OPERATOR_EXT_PBC_H

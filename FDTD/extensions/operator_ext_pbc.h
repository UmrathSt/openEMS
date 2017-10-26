#ifndef OPERATOR_EXT_PBC_H
#define OPERATOR_EXT_PBC_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Operator_Ext_Pbc : public Operator_Extension
{
    friend class Engine_Ext_Pbc;
    friend class Operator_Extension;
    friend class Operator;
public:
    Operator_Ext_Pbc(Operator* op);
    ~Operator_Ext_Pbc();

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

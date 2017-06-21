#ifndef OPERATOR_EXT_PBC_H
#define OPERATOR_EXT_PBC_H

#include "operator_extension.h"
#include "FDTD/operator.h"

class Pbc;

class Operator_Ext_Pbc : public Operator_Extension
{
    friend class Engine_Ext_Pbc;
    friend class Operator;
public:
    Operator_Ext_Pbc(Operator* op);
    ~Operator_Ext_Pbc();

    virtual Operator_Extension* Clone(Operator* op);
    // sets the phase difference between opposite sides of the periodic structure
    // Example: periodicity (exp(i * k * r) in x-direction only -> kparallel = {1, 0, 0};
    void SetKParallel(double kparallel[3]);

    virtual bool BuildExtension();

    virtual Engine_Extension* CreateEngineExtention();

    virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
    virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}
    virtual bool IsMPISave() const {return true;}

    virtual string GetExtensionName() const {return string("PeriodicBoundaryCondition Extension");}

    virtual void ShowStat(ostream &ostr) cons

protected:
    Operator_Ext_Pbc(Operator* op, Operator_Ext_Pbc* op_ext);
};


#endif // OPERATOR_EXT_PBC_H

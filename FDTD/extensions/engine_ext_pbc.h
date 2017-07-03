#ifndef ENGINE_EXT_PBC_H
#define ENGINE_EXT_PBC_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

class Operator_Ext_Pbc;

class Engine_Ext_Pbc : public Engine_Extension
{
public:
    Engine_Ext_Pbc(Operator_Ext_Pbc* op_ext);
    virtual ~Engine_Ext_Pbc();

    virtual void Apply2Voltages();
    virtual void Apply2Current();

protected:
    Operator_Ext_Pbc* m_Op_Pbc;
};


#endif // ENGINE_EXT_PBC_H

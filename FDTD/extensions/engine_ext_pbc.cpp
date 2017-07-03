/* engine extension for periodic boundary conditions
*/
#include "engine_ext_pbc.h"
#include "operator_extension.h"

#include "FDTD/engine.h"

Engine_Ext_Pbc::Engine_Ext_Pbc(Operator_Ext_Pbc* op_ext) : Engine_Extension(op_ext)
{
    m_Op_Pbc = op_ext;
    m_direction
}

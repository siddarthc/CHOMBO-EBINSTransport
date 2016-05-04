#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBScalarAdvectBCFactory.H"
#include "EBScalarAdvectBC.H"

#include "NamespaceHeader.H"

/******************/
EBScalarAdvectBCFactory::
EBScalarAdvectBCFactory(int                                    a_inflowDir,
                        Real                                   a_scalarInflowValue,
                        RefCountedPtr<PoiseuilleInflowBCValue> a_injectionBCFunc)
  :EBPhysIBCFactory()
{
  m_inflowDir         = a_inflowDir;
  m_scalarInflowValue = a_scalarInflowValue;
  m_injectionBCFunc   = a_injectionBCFunc;
}
/******************/
EBScalarAdvectBCFactory::
 ~EBScalarAdvectBCFactory()
{;}
/******************/
EBPhysIBC*
EBScalarAdvectBCFactory::
create() const
{
  EBScalarAdvectBC* retval = new EBScalarAdvectBC(m_inflowDir, m_scalarInflowValue, m_injectionBCFunc);

  return static_cast<EBPhysIBC*>(retval);
}
/******************/

#include "NamespaceFooter.H"

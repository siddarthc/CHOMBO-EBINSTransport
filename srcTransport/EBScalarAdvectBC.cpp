#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBScalarAdvectBC.H"
#include "EBLevelDataOps.H"
#include "EBArith.H"
#include "DirichletPoissonEBBC.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

extern Real g_simulationTime;

/****************************/
EBScalarAdvectBC::
~EBScalarAdvectBC()
{
  if (m_advVelPtr != NULL) delete m_advVelPtr;
}
/****************************/
EBScalarAdvectBC::
EBScalarAdvectBC(int a_inflowDir,
                 Real a_scalarInflowValue,
                 RefCountedPtr<PoiseuilleInflowBCValue> a_injectionBCFunc)
  :EBPhysIBC()
{
  m_inflowDir = a_inflowDir;
  m_scalarInflowValue = a_scalarInflowValue;
  m_injectionBCFunc = a_injectionBCFunc;
  m_isAdvVelSet = false;
}
/****************************/
void
EBScalarAdvectBC::
define(const ProblemDomain&  a_domain,
       const RealVect&       a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx;
  CH_assert(Abs(m_dx[0]-m_dx[1]) < 1.e-9);
  m_isDefined = true;
}
/****************************/
void
EBScalarAdvectBC::
fluxBC(EBFluxFAB&            a_facePrim,
       const EBCellFAB&      a_Wcenter,
       const EBCellFAB&      a_Wextrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_faceDir)
{
//  a_facePrim.setVal(1.0); 
  CH_assert(m_isDefined);
  CH_assert(!m_domain.isPeriodic(a_faceDir));
  CH_assert(m_injectionBCFunc != NULL);

  //gets face centered region of data
  Box FBox = a_facePrim[a_faceDir].getRegion();
  Box cellBox = a_faceBox;
  CH_assert(a_facePrim[a_faceDir].nComp() == 1);  

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_faceDir,isign);

  //figure out if we are on a wall and if its an inflow and the inflow value
  bool isInflow   = ((a_side == Side::Lo) && (a_faceDir==m_inflowDir));
  bool isOutflow  = ((a_side == Side::Hi) && (a_faceDir==m_inflowDir));

    // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box bndryBox = adjCellBox(cellBox, a_faceDir, a_side, 1);

      // Shift things to all line up correctly
      bndryBox.shift(a_faceDir,-isign);

      IntVectSet ivs(bndryBox);
      Real scalVal, velBC;
      for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_faceDir, a_side);
          for (int iface= 0; iface < bndryFaces.size(); iface++)
            {
              //this all will work for the case of spatially varying inflow
              //for inhomogeneous other faces, this is wrong.
              const FaceIndex& face = bndryFaces[iface];
              if (isOutflow) // assuming injection is from the outflow face (injectibe into counter flow)
                {
                  // figure out if you are inside or outside the tube
                  // inside the tube scalar injection val = 1.0; outside, extrap
 
                  RealVect prob_lo = RealVect::Zero;
                  const RealVect tubeCenter = m_injectionBCFunc->getTubeCenter();
                  Real tubeRadius = m_injectionBCFunc->getTubeRadius();
                  bool isInsideTube = false;
                  const RealVect loc  = EBArith::getFaceLocation(face, m_dx, prob_lo);
                  CH_assert(loc[m_inflowDir] == tubeCenter[m_inflowDir]);
                  Real radius = m_injectionBCFunc->getRadius(loc);
// start hardwire
                  Real wallThickness = 0.075;
                  ParmParse pp;
                  pp.get("wall_thickness",wallThickness);
                  Real bcFactor = 0.5; // % of outflow face to be set to inflow condition
                  tubeRadius += wallThickness * bcFactor /2.0;

                  if (radius <= tubeRadius)                 
                    {
                      isInsideTube = true;
                    }

                  isInsideTube ? scalVal = m_scalarInflowValue : scalVal = a_Wextrap(vof, 0);
                } // end if outflow
              else if (isInflow)
                {
                  scalVal = 0.0;
                }
              else  // solid wall
                {
                  scalVal = 0.0;  
                }
               a_facePrim[a_faceDir](face,0) = scalVal;
             } // end face loop
          }  // end vofit loop
       }          
}
/****************************/
void
EBScalarAdvectBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
  EBLevelDataOps::setToZero(a_conState);
}
/****************************/
void
EBScalarAdvectBC::
setAdvVel(const EBFluxFAB& a_advVel)
{
  m_isAdvVelSet = true;
  m_advVelPtr = &a_advVel;
}
/****************************/

#include "NamespaceFooter.H"

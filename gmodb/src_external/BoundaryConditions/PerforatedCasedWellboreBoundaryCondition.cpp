//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2014, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//
//  Randolph Settgast		Stuart Walsh
//  Scott Johnson		Pengcheng Fu
//  Joshua White
//
//  LLNL-CODE-656690
//  GMOD-B, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GMOD-B. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//  Please also read "Additional BSD Notice" below.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the 
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Additional BSD Notice
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file BoundaryConditions.cpp
 * @author walsh24
 * @date December 5, 2011
 */

#include "PerforatedCasedWellboreBoundaryCondition.h"
#include "Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "MPI_Communications/Communication.h"

#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

namespace {
  const realT pi = 3.14159265358979323846;

  realT CalculatePerforationAlpha(realT perforationDiameter,
		                          realT perforationChannelFriction,
		                          realT casingThickness,
		                          realT perforationK,
		                          realT rho,
		                          realT perforationCountPerCluster){

      const realT area = perforationDiameter * perforationDiameter * pi * 0.25;
      realT perforationAlpha = perforationChannelFriction * casingThickness;
      perforationAlpha /= 2 * perforationDiameter;
      perforationAlpha += perforationK;
      perforationAlpha *= rho;
      perforationAlpha /= area * area * perforationCountPerCluster;
      return perforationAlpha;
      //this is the 1/perm for the perf
  };

}



inline realT GetIorFirst(const rArray1d& array, unsigned i){
  return (array.size() == 1)? array[0]: array[i];
}

inline localIndex GetIorFirst(const lArray1d& array, unsigned i){
  return (array.size() == 1)? array[0]: array[i];
}

inline void Embiggen(localIndex& count, localIndex n){
  if(count > 1 && n > 1 && n!= count){
	  throw GPException(
	          "PerforatedCasedWellboreBoundaryCondition: Inconsistent perforation cluster data lengths: " + toString(count) +  " and " + toString(n) );

  }
  if(n>count) count = n;
}

////////////////////////////
//
/// Perforated wellbore boundary condition
/**
 * @author johnson346
 * @brief Simple wellbore hydraulics model
 *
 **/

PerforatedCasedWellboreBoundaryCondition::PerforatedCasedWellboreBoundaryCondition(
    TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm) :
    RadialHydraulicPressureBoundaryCondition(hdn, pm, false),
    m_setNamesHydro(),
    m_radialDistance(-1.0),
    m_mass(0.0),
    m_volume(-1.0),
    m_K0(0.0),
    m_rho0(0.0),
    m_pressureCap(0.0),
    m_pipeFrictionFactor(0.0),
    m_perforationAlpha(0.0),
    m_perforationSnapDistance(0.0),
    m_locationsSet(false),
    m_perforationLocations(),
    m_initialGuess()
#if GPAC_MPI
    ,
    m_globalToNeighborRank(),
    m_globalToLocal(),
    m_sizes()
#endif
{
  m_fieldName = "perforatedCasedWellbore";

  m_axialDistance *= 2.0;
  m_radialDistance = hdn->GetAttributeOrDefault<realT>("radius", -1.0);
  m_K0 = hdn->GetAttributeOrDefault<realT>("bulkModulus", 2e9);
  m_rho0 = hdn->GetAttributeOrDefault<realT>("referenceDensity", 1000);
  m_pressureCap = hdn->GetAttributeOrDefault<realT>("pressureCap",
                                                    std::numeric_limits<realT>::max());
  if (m_pressureCap > m_K0)
    throw GPException(
        "PerforatedCasedWellboreBoundaryCondition: Pressure cap should be less than bulk modulus");

  {
    //PERFORATION CALCS
	  /*
    const realT perforationChannelFriction = hdn->GetAttributeOrDefault<realT>(
        "perforationChannelFriction", 0.6);
    const realT casingThickness = hdn->GetAttributeOrDefault<realT>("casingThickness", 0.6);
    const realT perforationK = hdn->GetAttributeOrDefault<realT>("perforationK", 60.0);
    const realT perforationDiameter = hdn->GetAttributeOrDefault<realT>("perforationDiameter",
                                                                        0.03);
    const localIndex perforationCountPerCluster = hdn->GetAttributeOrDefault<localIndex>(
        "perforationCountPerCluster", 10);
    */

	    rArray1d perforationChannelFrictions;
	    rArray1d casingThicknesses;
	    rArray1d perforationKs;
	    rArray1d perforationDiameters;
	    lArray1d perforationCountsPerCluster;

	    perforationChannelFrictions = hdn->GetAttributeVectorOrDefault<realT>(
	    	        "perforationChannelFriction", std::vector<realT>(1,0.6));
	    casingThicknesses = hdn->GetAttributeVectorOrDefault<realT>("casingThickness", std::vector<realT>(1,0.6));
	    perforationKs = hdn->GetAttributeVectorOrDefault<realT>("perforationK", std::vector<realT>(1,60.0));
	    perforationDiameters = hdn->GetAttributeVectorOrDefault<realT>("perforationDiameter", std::vector<realT>(1,0.03));
	    perforationCountsPerCluster = hdn->GetAttributeVectorOrDefault<localIndex>(
	        "perforationCountPerCluster", std::vector<localIndex>(1,10));


    /*
    const realT area = perforationDiameter * perforationDiameter * pi * 0.25;
    m_perforationAlpha = perforationChannelFriction * casingThickness;
    if(isZero(perforationDiameter))
      throw GPException("Impermeable perforations specified: perforation diameter zero!");
    m_perforationAlpha /= 2 * perforationDiameter;
    m_perforationAlpha += perforationK;
    m_perforationAlpha *= m_rho0;
    if(perforationCountPerCluster == 0)
      throw GPException("Impermeable perforations specified: zero perforations per cluster!");
    m_perforationAlpha /= area * area * perforationCountPerCluster;
    */
    /*
    m_perforationAlpha = CalculatePerforationAlpha(perforationDiameter,
                                                   perforationChannelFriction,
                                                   casingThickness,
                                                   perforationK,
                                                   m_rho0,
                                                   perforationCountPerCluster);
    */

    localIndex numberOfPerforationClusters = 1;
    Embiggen(numberOfPerforationClusters,perforationDiameters.size());
    Embiggen(numberOfPerforationClusters,perforationChannelFrictions.size());
    Embiggen(numberOfPerforationClusters,casingThicknesses.size());
    Embiggen(numberOfPerforationClusters,perforationKs.size());
    Embiggen(numberOfPerforationClusters,perforationCountsPerCluster.size());

    for(localIndex i =0; i < numberOfPerforationClusters; ++i){
      realT perforationDiameter  = GetIorFirst(perforationDiameters,i);
      realT perforationChannelFriction = GetIorFirst(perforationChannelFrictions,i);
      realT casingThickness = GetIorFirst(casingThicknesses,i);
      realT perforationK     = GetIorFirst(perforationKs,i);
      localIndex perforationCountPerCluster = GetIorFirst(perforationCountsPerCluster,i);




      if(perforationCountPerCluster == 0)
        throw GPException("Impermeable perforations specified: zero perforations per cluster!");
      if(isZero(perforationDiameter))
        throw GPException("Impermeable perforations specified: perforation diameter zero!");


      realT alpha = CalculatePerforationAlpha(perforationDiameter,
                    perforationChannelFriction,
                    casingThickness,
                    perforationK,
                    m_rho0,
                    perforationCountPerCluster);

      if(alpha <= 0)
        throw GPException("Frictionless perfs specified!");

      m_perforationAlphas.push_back(alpha);
    }

    m_perforationAlpha = m_perforationAlphas[0]; // for backwards compatability for the time being
    //this is the 1/perm for the perf

    //PIPE CALCS
    const realT pipeFriction = hdn->GetAttributeOrDefault<realT>("pipeFriction", 0.6);
    if(pipeFriction <= 0)
      throw GPException("Frictionless pipe specified: pipe friction is zero or negative!");
    if(m_radialDistance <= 0)
      throw GPException("Impermeable pipe specified: pipe diameter is zero or negative!");
    const realT areaPipe = m_radialDistance * m_radialDistance * pi;
    m_pipeFrictionFactor = pipeFriction;
    m_pipeFrictionFactor /= 4 * m_radialDistance;
    m_pipeFrictionFactor *= m_rho0;
    m_pipeFrictionFactor /= areaPipe * areaPipe;

    //m_pipeFrictionFactor * L is the "pipeAlpha" or 1/perm
  }
  const realT lEff = hdn->GetAttributeOrDefault<realT>("effectiveBoreholeLength", m_axialDistance);

  //set the distance between perfs, beneath which all perforations are treated the same
  m_perforationSnapDistance = hdn->GetAttributeOrDefault<realT>("perforationSnap", 0.3);

  m_volume = lEff > 0 ? lEff : 0;
  m_volume *= m_radialDistance > 0 ? (m_radialDistance * m_radialDistance) : 0;
  m_volume *= pi;

  m_volume = hdn->GetAttributeOrDefault<realT>("volume", m_volume);

  m_mass = hdn->GetAttributeOrDefault<realT>("initialMass", m_rho0 * m_volume);

  m_setNamesHydro = hdn->GetStringVector("setnamesHydro");
}

// Individual solvers should use this function rather than calculate the traction directly in the solver.
// This allows the traction bc to be developed separately from the individual solvers
R1Tensor PerforatedCasedWellboreBoundaryCondition::GetTractionOnFace(PhysicalDomainT& domain,
                                                                     const lSet::const_iterator& fc,
                                                                     realT& time)
{
  R1Tensor traction(0);
  //FOR NOW, LET'S ASSUME A CASED WELL-BORE WHERE THE CASING SUPPORTS THE LOAD OF THE FLUID PRESSURE
#if 0
  if(m_axialDistance <= 0 || m_volume <= 0)
  return traction;

  R1Tensor position;
  {
    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *fc, position );
    position -= m_origin;

    const realT distance = Dot(position, m_axis);
    if(fabs(distance) > m_axialDistance)
    return traction;

    R1Tensor n(m_axis);
    n *= distance;
    position -= n;
  }
  //We can assume the face set has been declared to describe the well bore wall
  //  const realT radius = position.L2_Norm();
  //  if(radius > m_radialDistance)
  //    return traction;

  traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
  traction *= -GetPressure();
#endif
  return traction;
}

void PerforatedCasedWellboreBoundaryCondition::Apply(PhysicalDomainT& domain, realT time, realT dt)
{
  //Do the initialization if necessary
  //the 2nd argument will force setting if "true"
  SetPerforationLocations(domain, false);

  //Buffer the pressures
  rArray1d fracturePressures, fractureFluxes, facePerms;
  FaceManagerT& fm = domain.m_feFaceManager;
  int rank = ExtractFaceData(fm, fracturePressures, facePerms);

#if 0
  //DEBUG
  std::cout << "# of perfs: " << m_perforationLocations.size() << std::endl;
  for(Array1dT<ValueGUID>::const_iterator it = m_perforationLocations.begin(); it != m_perforationLocations.end(); ++it)
  {
    std::cout << "  perf location (" << it->m_index << ") = " << it->m_x << std::endl;
  }
  std::cout << "# of pressures: " << fracturePressures.size() << std::endl;
  for(rArray1d::const_iterator it = fracturePressures.begin(); it != fracturePressures.end(); ++it)
  {
    std::cout << "  pressure = " << *it << std::endl;
  }
  throw GPException("DEBUG");
#endif

  if (rank == 0)
    Solve(time, fracturePressures, facePerms, fractureFluxes);

#if GPAC_MPI
  //Fill data structures and return to owners
  FractureFluxesToAll(fm, dt, fractureFluxes);
#else
  {
    const rArray1d& fluidVolume  = fm.GetFieldData<FieldInfo::volume>();
    rArray1d& faceFluidMass = fm.GetFieldData<FieldInfo::mass>();
    rArray1d& faceFluidDensity = fm.GetFieldData<FieldInfo::density>();
    const std::map<globalIndex, localIndex> gtol = fm.m_globalToLocalMap;
    rArray1d::const_iterator itf = fractureFluxes.begin();
    for(Array1dT<ValueGUID>::iterator it = m_perforationLocations.begin(); it != m_perforationLocations.end(); ++it, ++itf)
    {
      std::map<globalIndex, localIndex>::const_iterator it2 = gtol.find(it->m_index);
      if(it2 == gtol.end())
        throw GPException("Cannot find request entry!!");
      const localIndex a = it2->second;
      faceFluidMass[a] += (*itf) * dt * m_rho0;
      faceFluidDensity[a] = faceFluidMass[a]/fluidVolume[a];
    }
  }
}
#endif
}

int PerforatedCasedWellboreBoundaryCondition::Fill(const rArray1d& fracturePressures,const rArray1d& facePerms,
                                                   rArray1d& betas, rArray1d& xs, rArray1d& ps, rArray1d& fp,
                                                   iArray1d& ns)
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //We need to reduce the problem by agglomerating all fracture faces into clusters that are collocated
  //(or close to it) along the wellbore:
  //
  //This is necessary, as it both reduces the degrees-of-freedom and removes infinite permeability terms.
  //
  //The result is a list of the beta terms (betas),
  //                        the agglomerated cluster positions (xs),
  //                        the agglomerated cluster fracture pressures (ps), and
  //                        the number of fracture faces per agglomerated cluster (ns)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  realT lastPosition = 0; //m_perforationLocations.begin()->m_x;
  int nlocs = 1;

  //for each pressure and perf location, check to see whether the cluster is new or pre-existing
  //if pre-existing add the pressure to the sum (average later)

  rArray1d::const_iterator it1 = fracturePressures.begin();
  rArray1d::const_iterator it2 = facePerms.begin();
  for (Array1dT<ValueGUID>::const_iterator it = m_perforationLocations.begin();
      it != m_perforationLocations.end(); ++it, ++it1, ++it2)
  {
    //get the distance along between the perf and the last (or the origin if it is the first)  
    const realT dx = it->m_x - lastPosition;
    if (dx > m_perforationSnapDistance)
    {
      if (nlocs > 1)
      {
        //if you need to unwind the last agglomeration, do it here
        ns.back() = nlocs;
        ps.back() /= ns.back();
        fp.back() /= ns.back();
      }
      nlocs = 1;

      //you do not need the first length OR beta if this is a flux BC
      if (it != m_perforationLocations.begin())
      {
        xs.push_back(dx);
        const realT beta = m_pipeFrictionFactor * dx;
        betas.push_back(beta);
      }
      ps.push_back(*it1);
      fp.push_back(*it2);
      ns.push_back(nlocs);
    }
    else
    {
      if (it == m_perforationLocations.begin())
        throw GPException(
            "Cannot have the first perforation location be at or near the injection point!!");
      ps.back() += *it1;
      fp.back() += *it2;
      ++nlocs;
    }
    lastPosition = it->m_x;
  }

  if (nlocs > 1)
  {
    //if you need to unwind the last agglomeration, do it here
    ns.back() = nlocs;
    ps.back() /= ns.back();
    fp.back() /= fp.back();
  }
  return 0;
}

void PerforatedCasedWellboreBoundaryCondition::Solve(const realT time,
                                                     const rArray1d& fracturePressures,
                                                     const rArray1d& facePerms,
                                                     rArray1d& fractureFluxes)
{
	//////////////////////////////////////////////////////////////
  //Get the wellbore flow rate: put into the m_value variable
  if (!(m_timeTableName.empty()))
  {
    rArray1d t(1);
    t[0] = time;
    const realT tableval = TableManager::Instance().LookupTable<1>(m_timeTableName, t);
    m_value = m_scale * tableval;
  }
  else if (!(m_functionName.empty()))
  {
    m_value = m_scale * (*m_function)(time);
  }

  //////////////////////////////////////////////////////////////
  //Calculate the agglomerated cluster parameters (betas, xs, ps, ns)
  rArray1d betas, xs, ps, fp;
  iArray1d ns;
  Fill(fracturePressures, facePerms, betas, xs, ps, fp, ns);
#if 0
  //DEBUG
  std::cout << "# of locations: " << xs.size() << std::endl;
  for(rArray1d::const_iterator it = xs.begin(); it != xs.end(); ++it)
  {
    std::cout << "  x = " << *it << std::endl;
  }
  std::cout << "# of betas: " << betas.size() << std::endl;
  for(rArray1d::const_iterator it = betas.begin(); it != betas.end(); ++it)
  {
    std::cout << "  beta = " << *it << std::endl;
  }
  std::cout << "# of pressures: " << ps.size() << std::endl;
  for(rArray1d::const_iterator it = ps.begin(); it != ps.end(); ++it)
  {
    std::cout << "  p = " << *it << std::endl;
  }
  std::cout << "# of fracs per location entries: " << ns.size() << std::endl;
  for(iArray1d::const_iterator it = ns.begin(); it != ns.end(); ++it)
  {
    std::cout << "  f per loc = " << *it << std::endl;
  }
  //throw GPException("DEBUG");
#endif

  //////////////////////////////////////////////////////////////
  //Solve the agglomerated system

  //Call the TrilinosUtilities non-linear solve
  NonlinearSolve(ps, fp, m_perforationAlpha, betas, m_value);

  //Set the new flux values
  {
    fractureFluxes.resize(fracturePressures.size());

    rArray1d::const_iterator itqf = m_initialGuess.begin();
    iArray1d::const_iterator itns = ns.begin();
    realT qcurr = (*itqf)/(*itns);
    int iii = 1;
    for (rArray1d::iterator it = fractureFluxes.begin(); it != fractureFluxes.end(); ++it)
    {
      //put the same flux into any perf region that is within this cluster
      *it = qcurr;

      //advance the cursor only if this flux is part of a new cluster
      if (iii >= *itns)
      {
        iii = 1;
        ++itqf;
        ++itns;
	if(itns != ns.end())
	  qcurr = (*itqf)/(*itns);
      }
      else
      {
        ++iii;
      }
    }
  }
}

void PerforatedCasedWellboreBoundaryCondition::NonlinearSolve(const rArray1d& fracturePressures, const rArray1d& facePerms, realT alpha,
                                                              const rArray1d& beta, realT qwellbore,
                                                              const bool printFlag)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

#ifdef GPAC_MPI
  Epetra_MpiComm CommWorld(MPI_COMM_WORLD);
#else
  Epetra_SerialComm CommWorld;
#endif

  // The example problem is so small that we must run it on only one
  // process.  However, people might run this example code in MPI with
  // any number of processes.  We handle this by using a communicator
  // containing only one MPI process, and quieting all processes but
  // Proc 0 (with respect to MPI_COMM_WORLD).
  if (CommWorld.MyPID() == 0)
  {
    // Set up the problem interface.  Your application will define
    // its own problem interface.  SimpleProblemInterface is our
    // example interface, which you can use as a model.
    //
    // Our SimpleProblemInterface makes a deep copy of the initial
    // guess and exact solution vectors.
    RCP<PerforatedCasedWellboreProblemInterfaceQD> interface = rcp(
        new PerforatedCasedWellboreProblemInterfaceQD());


    interface.get()->Initialize(fracturePressures, facePerms, m_perforationAlphas, beta, m_value);

#ifdef GPAC_MPI
    Epetra_MpiComm Comm (MPI_COMM_SELF);
#else
    Epetra_SerialComm Comm;
#endif

    const bool getInitialGuess = m_initialGuess.size() != (int)interface.get()->Size();
    if(getInitialGuess)
      m_initialGuess.resize(interface.get()->Size(),0.0);
    Epetra_Map Map((int)m_initialGuess.size(), 0, Comm);
    Epetra_Vector xx(Map);

    //set the inital guess either from the last solution
    //or from the estimator
    if(getInitialGuess)
    {
      localIndex i = 0;
      interface.get()->InitialGuess(xx);
      for(rArray1d::iterator it = m_initialGuess.begin(); it != m_initialGuess.end(); ++it, ++i)
        *it = xx[i];
    }
    else
    {
      localIndex i = 0;
      for(rArray1d::const_iterator it = m_initialGuess.begin(); it != m_initialGuess.end(); ++it, ++i)
        xx[i] = *it;
    }

    // Create the top-level parameter list to control NOX.
    //
    // "parameterList" (lowercase initial "p") is a "nonmember
    // constructor" that returns an RCP<ParameterList> with the
    // given name.
    RCP < ParameterList > params = parameterList("NOX");

    // Tell the nonlinear solver to use line search.
    params->set("Nonlinear Solver", "Line Search Based");

    //
    // Set the printing parameters in the "Printing" sublist.
    //
    ParameterList& printParams = params->sublist("Printing");
    printParams.set("MyPID", Comm.MyPID());
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);

    // Set verbose=true to see a whole lot of intermediate status
    // output, during both linear and nonlinear iterations.
    const bool verbose = printFlag;
    if (verbose)
    {
      printParams.set(
          "Output Information",
          NOX::Utils::OuterIteration + NOX::Utils::OuterIterationStatusTest + NOX::Utils::InnerIteration + NOX::Utils::Parameters + NOX::Utils::Details + NOX::Utils::Warning);
    }
    else
    {
      printParams.set("Output Information", NOX::Utils::Warning);
    }

    //
    // Set the nonlinear solver parameters.
    //

    // Line search parameters.
    ParameterList& searchParams = params->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Parameters for picking the search direction.
    ParameterList& dirParams = params->sublist("Direction");
    // Use Newton's method to pick the search direction.
    dirParams.set("Method", "Newton");

    // Parameters for Newton's method.
    ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    //
    // Newton's method invokes a linear solver repeatedly.
    // Set the parameters for the linear solver.
    //
    ParameterList& lsParams = newtonParams.sublist("Linear Solver");

    // Use Aztec's implementation of GMRES, with at most 800
    // iterations, a residual tolerance of 1.0e-4, with output every
    // 50 iterations, and Aztec's native ILU preconditioner.
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", 800);
    lsParams.set("Tolerance", 1e-8);
    lsParams.set("Output Frequency", 50);
    lsParams.set("Aztec Preconditioner", "ilu");

    //
    // Build the Jacobian matrix.
    //
    RCP < Epetra_CrsMatrix > A = rcp(new Epetra_CrsMatrix(Copy, Map, (int)m_initialGuess.size()));
    interface.get()->constructJacobian(xx, *A.get());

    // Our SimpleProblemInterface implements both Required and
    // Jacobian, so we can use the same object for each.
    RCP < NOX::Epetra::Interface::Required > iReq = interface;
    RCP < NOX::Epetra::Interface::Jacobian > iJac = interface;

    RCP < NOX::Epetra::LinearSystemAztecOO > linSys = rcp(
        new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iJac, A, xx));

    // Need a NOX::Epetra::Vector for constructor.
    NOX::Epetra::Vector noxInitGuess(xx, NOX::DeepCopy);
    RCP < NOX::Epetra::Group > group = rcp(
        new NOX::Epetra::Group(printParams, iReq, noxInitGuess, linSys));

    //
    // Set up NOX's iteration stopping criteria ("status tests").
    //
    // ||F(X)||_2 / N < 1.0e-4, where N is the length of F(X).
    //
    // NormF has many options for setting up absolute vs. relative
    // (scaled by the norm of the initial guess) tolerances, scaling
    // or not scaling by the length of F(X), and choosing a
    // different norm (we use the 2-norm here).
    RCP < NOX::StatusTest::NormF > testNormF = rcp(new NOX::StatusTest::NormF(1.0e-8));

    // At most 20 (nonlinear) iterations.
    RCP < NOX::StatusTest::MaxIters > testMaxIters = rcp(new NOX::StatusTest::MaxIters(20));

    // Combine the above two stopping criteria (normwise
    // convergence, and maximum number of nonlinear iterations).
    // The result tells NOX to stop if at least one of them is
    // satisfied.
    RCP < NOX::StatusTest::Combo > combo = rcp(
        new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters));

    // Create the NOX nonlinear solver.
    RCP < NOX::Solver::Generic > solver = NOX::Solver::buildSolver(group, combo, params);

    // Solve the nonlinear system.
#if 0
    NOX::StatusTest::StatusType status = solver->solve();
#else
    solver->solve();
#endif

    // Print the result.
    //
    // For this particular example, Comm contains only one MPI
    // process.  However, we check for Comm.MyPID() == 0 here just
    // so that the example is fully general.  (If you're solving a
    // larger nonlinear problem, you could safely use the code
    // below.)
    if (Comm.MyPID() == 0)
    {
      cout << endl << "-- Parameter List From Solver --" << endl;
      solver->getList().print(cout);
    }

    // Get the Epetra_Vector with the final solution from the solver.
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());

    xx = dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX()).getEpetraVector();

    if (Comm.MyPID() == 0)
    {
      cout << "Computed solution: " << endl;
    }
    // Epetra objects know how to print themselves politely when
    // their operator<<(std::ostream&) is invoked on all MPI
    // process(es) in the communicator to which they are associated.
    cout << xx;
    {
      localIndex i = 0;
      for(rArray1d::iterator it = m_initialGuess.begin(); it != m_initialGuess.end(); ++it, ++i)
        *it = xx[i];
    }
  }
}

//void PerforatedCasedWellboreBoundaryCondition::LinSolve(const realT time,
//                                                        const rArray1d& fracturePressures,
//                                                        rArray1d& fractureFluxes)
//{
//  //CHECK OUT NOX PACKAGE
//  //didasko/examples/nox/ex1.cpp
//  //Allocate and initialize
//  const localIndex size = 3 * m_perforationLocations.size() - 1;
//  rArray2d A;
//  A.resize2(size, size);
//  rArray1d x(size), b(size);
//  b = 0;
//  x = 0;
//  A = 0;
//
//  //Fill
//  //Solve the system: assume incompressible flow in a simple well-bore, where perforations are clustered
//  //
//  //The system is 3*N-1 unknowns (N=number of perforation clusters): IN ORDER
//  //  N pressures in the well-bore at the points where the clusters are placed
//  //  N fluxes through the casing into the initiated fractures
//  //  N-1 fluxes along the well-bore between clusters
//  //
//  //The equations can be formulated as:
//  // (1: i=1..N)     q_w,i + q_f,i - q_w,i-1           = 0
//  // (2: i=1..N)     p_w,i - alpha_i * |q_f,i| * q_f,i = p_f,i
//  // (3: i=1..N-1)   alpha_i * |q_f,i| * q_f,i
//  //                 - beta_i * |q_w,i| * q_w,i
//  //                 - alpha_i+1 * |q_f,i+1| * q_f,i+1 = p_f,i+1 - p_f,i
//
//  //FILL b
//  for (localIndex i1 = 0; i1 < size; i1++)
//  {
//    const localIndex i2 = i1 + size;
//    b[i2] = fracturePressures[i1];
//    if (i1 > 0)
//    {
//      const localIndex i3 = i2 + size - 1;
//      b[i3] = fracturePressures[i1] - fracturePressures[i1 - 1];
//    }
//    else
//    {
//      if (!(m_timeTableName.empty()))
//      {
//        rArray1d t(1);
//        t[0] = time;
//        const realT tableval = TableManager::Instance().LookupTable<1>(m_timeTableName, t);
//        m_value = m_scale * tableval;
//      }
//      else if (!(m_functionName.empty()))
//      {
//        m_value = m_scale * (*m_function)(time);
//      }
//      b[i1] = m_value;
//    }
//  }
//
//  //FILL A
//  for (localIndex i1 = 0; i1 < size; i1++)
//  {
//    const localIndex iqf = i1;
//    const localIndex iqw = i1 + size;
//    const localIndex ipw = iqw + size - 1;
//
//    // (1: i=1..N)     q_w,i + q_f,i - q_w,i-1           = 0
//    A(i1, iqw) = 1.0;
//    A(i1, iqf) = 1.0;
//    A(i1, iqw) = 1.0;
//    if (i1 > 0)
//      A(i1, iqw - 1) = -1.0;
//
//    // (2: i=1..N)     p_w,i - alpha_i * |q_f,i| * q_f,i = p_f,i
//    const localIndex i2 = i1 + size;
//    A(i2, ipw) = 1.0;
//    A(i2, iqf) = -m_perforationAlpha; // * fluxEstimate[i1];
//    if (i1 > 0)
//    {
//      const localIndex i3 = i2 + size - 1;
//      // (3: i=1..N-1)   alpha_i * |q_f,i| * q_f,i
//      //                 - beta_i * |q_w,i| * q_w,i
//      //                 - alpha_i+1 * |q_f,i+1| * q_f,i+1 = p_f,i+1 - p_f,i
//      A(i3, iqf - 1) = -m_perforationAlpha; // * fluxEstimate[i1-1];
//      A(i3, iqw - 1) = -m_pipeFrictionFactor * (m_perforationLocations[i1].m_x - m_perforationLocations[i1 - 1].m_x); // * fluxEstimate[iqw-1];
//      A(i3, iqf) = -m_perforationAlpha; // * fluxEstimate[i1];
//    }
//  }
//
//  //Solve
//  LinSolve_Local(A, x, b);
//}

int PerforatedCasedWellboreBoundaryCondition::ExtractFaceData(const FaceManagerT& fm,
                                                              rArray1d& pressures,
                                                              rArray1d& fperms )
{
  int rank = 0;
  pressures.clear();
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  pressures.resize(m_perforationLocations.size());
  fperms.resize(m_perforationLocations.size());
  FracturePressuresToRoot(fm, pressures, fperms);
#else
  {
    const rArray1d& faceFluidPressure = fm.GetFieldData<FieldInfo::pressure>();
    const rArray1d& facePerm = fm.GetFieldData<realT>("permeability");

    const std::map<globalIndex, localIndex> gtol = fm.m_globalToLocalMap();
    for(Array1dT<ValueGUID>::iterator it = m_perforationLocations.begin(); it != m_perforationLocations.end(); ++it)
    {
      std::map<globalIndex, localIndex>::const_iterator it2 = gtol.find(it->m_index);
      if(it2 == gtol.end())
        throw GPException("Cannot find request entry!!");
      pressures.push_back(faceFluidPressure[it2->second]);
      fperms.push_back(facePerm[it2->second]);
    }
  }
#endif
  return rank;
}

void PerforatedCasedWellboreBoundaryCondition::FracturePressuresToRoot(const FaceManagerT& fm,
                                                                       rArray1d& pressures,
                                                                       rArray1d& facePerms )
{
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  CommRegistry::commID commID = CommRegistry::genericComm01;

  const rArray1d& faceFluidPressure = fm.GetFieldData<FieldInfo::pressure>();
  const rArray1d& facePerm = fm.GetFieldData<realT>("permeability");

  const std::map<globalIndex, localIndex>& gtol = fm.m_globalToLocalMap;

  if (rank == 0)
  {
    //communicate if necessary
    if(size > 1)
    {
      Array1dT<Array1dT<ValueGUID> > receiveBuffers(size - 1);
      Array1dT<MPI_Request> inReq(size - 1);
      Array1dT<MPI_Status> inStat(size - 1);
      for (int i = 1; i < size; i++)
      {
        //COMM_C: RECEIVE FLUID PRESSURES
        const int rtag = CommRegistry::CommTag(i, 0, commID);
        const int rsize = m_sizes[i] * sizeof(ValueGUID);
        receiveBuffers[i-1].resize(m_sizes[i]);
        MPI_Irecv(receiveBuffers[i - 1].data(), rsize, MPI_CHAR, i, rtag, MPI_COMM_WORLD,
                  &inReq[i - 1]);
      }

      for (unsigned int count = 0; count < inReq.size(); ++count)
      {
        int neighborIndex;
        MPI_Waitany(inReq.size(), inReq.data(), &neighborIndex, inStat.data());

        //deal with received buffers
        Array1dT<ValueGUID>& rbuf = receiveBuffers[neighborIndex];
        int i = 0;
        for (Array1dT<ValueGUID>::iterator it = rbuf.begin(); it != rbuf.end(); ++it, ++i)
        {
          pressures[m_globalToLocal[it->m_index]] = it->m_x;
          facePerms[m_globalToLocal[it->m_index]] = it->m_y;
        }
      }
    }

    //set my own pressures
    for (Array1dT<ValueGUID>::iterator it = m_perforationLocations.begin();
        it != m_perforationLocations.end(); ++it)
    {
      if (m_globalToNeighborRank[it->m_index] == 0)
      {
        std::map<globalIndex, localIndex>::const_iterator it2 = gtol.find(it->m_index);
        if (it2 == gtol.end())
          throw GPException("Cannot find requested entry!!");
        pressures[m_globalToLocal[it->m_index]] = faceFluidPressure[it2->second];
        facePerms[m_globalToLocal[it->m_index]] = facePerm[it2->second];
      }
    }
  }
  else
  {
    MPI_Request outReq;
    MPI_Status outStat;
    const int stag = CommRegistry::CommTag(rank, 0, commID);
    const int outSize = m_perforationLocations.size() * sizeof(ValueGUID);
    Array1dT<ValueGUID> pressureBuffer;
    for (Array1dT<ValueGUID>::iterator it = m_perforationLocations.begin();
        it != m_perforationLocations.end(); ++it)
    {
      std::map<globalIndex, localIndex>::const_iterator it2 = gtol.find(it->m_index);
      if (it2 == gtol.end())
        throw GPException("Cannot find requested entry!!");
      pressureBuffer.push_back(ValueGUID(it->m_index, faceFluidPressure[it2->second], facePerm[it2->second] ) );
    }
    //COMM_C: SEND FLUID PRESSURES
    MPI_Isend(pressureBuffer.data(), outSize, MPI_CHAR, 0, stag, MPI_COMM_WORLD, &outReq);
    MPI_Wait(&outReq, &outStat);
  }
}

void PerforatedCasedWellboreBoundaryCondition::FractureFluxesToAll(FaceManagerT& fm, realT dt,
                                                                   const rArray1d& fluxes)
{
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  CommRegistry::commID commID = CommRegistry::genericComm01;

  const rArray1d& fluidVolume  = fm.GetFieldData<FieldInfo::volume>();
  rArray1d& faceFluidMass = fm.GetFieldData<FieldInfo::mass>();
  rArray1d& faceFluidDensity = fm.GetFieldData<FieldInfo::density>();

  const std::map<globalIndex, localIndex> gtol = fm.m_globalToLocalMap;

  if (rank == 0)
  {
    Array1dT<Array1dT<ValueGUID> > sendBuffers;
    Array1dT<MPI_Request> outReq;
    Array1dT<MPI_Status> outStat;
    if(size > 1)
    {
      sendBuffers.resize(size - 1);
      outReq.resize(size - 1);
      outStat.resize(size - 1);
      for (int i = 1; i < size; i++)
      {
        const int im1 = i - 1;
        for (std::map<globalIndex, localIndex>::const_iterator it = m_globalToNeighborRank.begin();
            it != m_globalToNeighborRank.end(); ++it)
          if (it->second == (localIndex)i)
            sendBuffers[im1].push_back(ValueGUID(it->first, fluxes[m_globalToLocal[it->first]], fluxes[m_globalToLocal[it->first]]));
        const int stag = CommRegistry::CommTag(i, 0, commID);
        const int ssize = m_sizes[i] * sizeof(ValueGUID);
        MPI_Isend(sendBuffers[im1].data(), ssize, MPI_CHAR, i, stag, MPI_COMM_WORLD,
                  &outReq[im1]);
      }
    }

    //set my own fluxes
    rArray1d::const_iterator itf = fluxes.begin();
    for (Array1dT<ValueGUID>::iterator it = m_perforationLocations.begin();
        it != m_perforationLocations.end(); ++it, ++itf)
    {
      if (m_globalToNeighborRank[it->m_index] == 0)
      {
        std::map<globalIndex, localIndex>::const_iterator it2 = gtol.find(it->m_index);
        if (it2 == gtol.end())
          throw GPException("Cannot find requested entry!!");
        const localIndex a = it2->second;
        faceFluidMass[a] += (*itf) * dt * m_rho0;
        faceFluidDensity[a] = faceFluidMass[a] / fluidVolume[a];
#if 1
        //DEBUG
        std::cout << "GLOBALINDEX " << it2->first << " FLUX " << *itf << std::endl;
#endif
      }
    }

    //wait for all to finish before releasing memory
    MPI_Waitall(outReq.size(), outReq.data(), outStat.data());
  }
  else
  {
    MPI_Request inReq;
    MPI_Status inStat;
    const int rtag = CommRegistry::CommTag(rank, 0, commID);
    const int rsize = m_perforationLocations.size() * sizeof(ValueGUID);
    Array1dT<ValueGUID> fluxBuffer(m_perforationLocations.size());
    MPI_Irecv(fluxBuffer.data(), rsize, MPI_CHAR, 0, rtag, MPI_COMM_WORLD, &inReq);
    MPI_Wait(&inReq, &inStat);

    //process my fluxes
    for (Array1dT<ValueGUID>::iterator it = fluxBuffer.begin(); it != fluxBuffer.end(); ++it)
    {
      std::map<globalIndex, localIndex>::const_iterator it2 = gtol.find(it->m_index);
      if (it2 == gtol.end())
        throw GPException("Cannot find request entry!!");
      const localIndex a = it2->second;
      faceFluidMass[a] += it->m_x * dt * m_rho0;
      faceFluidDensity[a] = faceFluidMass[a]  / fluidVolume[a];
#if 0
      //DEBUG
      std::cout << "GLOBALINDEX " << it2->first << " FLUX " << it->m_x << std::endl;
#endif
    }
  }
}


void PerforatedCasedWellboreBoundaryCondition::SetPerforationLocations(PhysicalDomainT& domain,
                                                                       bool force)
{
  if (!force && m_locationsSet)
    return;
  m_locationsSet = true;

#if GPAC_MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //(1) set my own m_perforationLocations (using ValueGUID) and sort, regardless of my rank
  const iArray1d& flowFaceTypes = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const iArray1d& ghostRanks = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const gArray1d ltog = domain.m_feFaceManager.m_localToGlobalMap;

#if GPAC_MPI
  m_globalToNeighborRank.clear();
  m_globalToLocal.clear();
  m_sizes.clear();
#endif
  m_perforationLocations.clear();

  //for each set composing the perforations
  for (sArray1d::const_iterator its = m_setNamesHydro.begin(); its != m_setNamesHydro.end();
      ++its)
  {
    //make sure the set exists
    std::map<std::string, lSet>::const_iterator setMap = domain.m_feFaceManager.m_Sets.find(*its);
    if (setMap == domain.m_feFaceManager.m_Sets.end())
      continue;

    //for each face in the set
    const lSet& faceSet = setMap->second;
    for (lSet::const_iterator fc = faceSet.begin(); fc != faceSet.end(); ++fc)
    {
      //if this is not a flow face or it is a ghost, skip it
      if(flowFaceTypes[*fc] < 0 || ghostRanks[*fc] >= 0)
        continue;

      //get the position and distance along the wellbore
      R1Tensor position;
      domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, *fc, position);
      position -= m_origin;
      const realT distance = Dot(position, m_axis);

      //if the perf is along the line segment of the wellbore, then insert it; otherwise, continue
      if (distance < 0 || distance > m_axialDistance)
        continue;

      //we now know that the face is in the set, owned by the process, is a flow face, and is along the wellbore segment
      //so we should be ready to add it to our list ... we just need to figure out where, since we will be sorting as we go
      bool ok = false;
      const globalIndex iglobal = ltog[*fc];
#if GPAC_MPI
      m_globalToNeighborRank[iglobal] = rank;
#endif
      for (Array1dT<ValueGUID>::iterator it2 = m_perforationLocations.begin();
          it2 != m_perforationLocations.end(); ++it2)
      {
        if (distance < (*it2).m_x)
        {
          m_perforationLocations.insert(
              it2, ValueGUID(iglobal, distance, distance));
          ok = true;
          break;
        }
      }
      if (!ok)
        m_perforationLocations.push_back(
            ValueGUID(iglobal, distance, distance));
    }//for each face in set
  }//for each set

#if GPAC_MPI

  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  if(size < 2)
//    return; //this step is unncessary

  CommRegistry::commID commID = CommRegistry::genericComm01;

  if (rank == 0)
  {
    //(2-root) initialize receive size buffer and gather sizes
    m_sizes.resize(size);
    {
      //COMM_A: GATHER BUFFER SIZES
      int sendArray = m_perforationLocations.size();
      MPI_Gather(&sendArray, 1, MPI_INT, m_sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    Array1dT<Array1dT<ValueGUID> > receiveBuffers(size - 1);
    Array1dT<MPI_Request> inReq(size - 1);
    Array1dT<MPI_Status> inStat(size - 1);
    for (int i = 1; i < size; i++)
    {
      //COMM_B: RECEIVE LOCATIONS
      const int im1 = i - 1;
      const int rtag = CommRegistry::CommTag(i, 0, commID);
      const int rsize = m_sizes[i] * sizeof(ValueGUID);
      receiveBuffers[im1].resize(m_sizes[i]);
      MPI_Irecv(receiveBuffers[im1].data(), rsize, MPI_CHAR, i, rtag, MPI_COMM_WORLD,
                &inReq[im1]);
    }

    //(3-root) set the map that tracks the receipt order for every received face
    //as the buffers roll in
    //(based on neighbor rank): m_globalToNeighborRank
    for (unsigned int count = 0; count < inReq.size(); ++count)
    {
      int neighborIndex;
      MPI_Waitany(inReq.size(), inReq.data(), &neighborIndex, inStat.data());
      //the returned index is the index in the inReq array ... this is rank-1
      ++neighborIndex;

      //deal with received buffers ... insert into positions as needed
      //might as do the insertion sort as buffers are rolling in and mind some
      //of the comm cost
      Array1dT<ValueGUID>& rbuf = receiveBuffers[neighborIndex - 1];
      int i = 0;
      for (Array1dT<ValueGUID>::iterator it = rbuf.begin(); it != rbuf.end(); ++it, ++i)
      {
        bool ok = false;
        m_globalToNeighborRank[it->m_index] = neighborIndex;
        for (Array1dT<ValueGUID>::iterator it2 = m_perforationLocations.begin();
            it2 != m_perforationLocations.end(); ++it2)
        {
          if ((*it).m_x < (*it2).m_x)
          {
            m_perforationLocations.insert(it2, *it);
            ok = true;
            break;
          }
        }
        if (!ok)
          m_perforationLocations.push_back(*it);
      }
    }

    //(4-root) set global to local map
    {
      localIndex i = 0;
      for (Array1dT<ValueGUID>::iterator it = m_perforationLocations.begin();
          it != m_perforationLocations.end(); ++it, ++i)
      {
        m_globalToLocal[it->m_index] = i;
      }
    }
  }
  else
  {
    //(2) initialize receive size buffer and send size
    {
      //COMM_A: GATHER BUFFER SIZES
      int rbuf;
      int sendarray = m_perforationLocations.size();
      MPI_Gather(&sendarray, 1, MPI_INT, &rbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    //(3) send my positions
    //COMM_B: SEND LOCATIONS
    MPI_Request outReq;
    MPI_Status outStat;
    const int stag = CommRegistry::CommTag(rank, 0, commID);
    const int outSize = m_perforationLocations.size() * sizeof(ValueGUID);
    MPI_Isend(m_perforationLocations.data(), outSize, MPI_CHAR, 0, stag, MPI_COMM_WORLD, &outReq);
    MPI_Wait(&outReq, &outStat);
  }
#endif
}

REGISTER_BoundaryCondition(PerforatedCasedWellboreBoundaryCondition)

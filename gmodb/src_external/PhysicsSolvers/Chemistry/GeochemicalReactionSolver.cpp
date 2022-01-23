//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)     Stuart Walsh(walsh24@llnl.gov)
//  Scott Johnson (johnson346@llnl.gov)        Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)           
//
//  All rights reserved.
//
//  This file is part of GPAC.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "GeochemicalReactionSolver.h"

#include "OneDReactionFront.h"

#include "ObjectManagers/UnitManager.h"
//#include "Utilities/ReducedRowEchelonForm.h"

GeochemicalReactionSolver::GeochemicalReactionSolver( const std::string& name,
                                                      ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

/*
 * 
 *   <GeochemicalReactionSolver name="grs" faceset="Zmax" temperature="200">    

        <SolidPhases names = "Quartz, Albite"
                     reactionfunctions="Quartz_ReactionRate, Albite_ReactionRate"
                     reactionfunctionvariables="H+"
                     equilibrium_constants = "Kq, Ka"    
                     />  

        <Species names = "H+, OH-, Ca+2, CaHCO3+, CaCO3, CO3-2, HCO3-, H2CO3, CO2" />

        <EquilibriumData
                 equilibrium_equations= 
                           "Kq = [SiO2];
                            Ka = [Al+3][Na+][SiO2]^3/[H+]^4"
                 logKdata=" Kq, -2.3823,  0.0,  0.0,  0.0,   0.0;
                            Ka, -5.95,    0.0,  0.0,  0.0,   0.0"
         />

      </GeochemicalReactionSolver>
 * 
 */
void GeochemicalReactionSolver::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{  
  SolverBase::ReadXML(hdn);

  int rank(0);
  #if GPAC_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  using namespace GPUnits;
  UnitManager& theUnitManager = UnitManager::Instance();
  
  // get species
  //////////////
  if (rank == 0)  std::cout << "Species" << std::endl;
  HierarchicalDataNode* speciesNode = hdn->GetChild("Species");
  {
    m_species.names = speciesNode->GetStringVector("names", ",");
    m_species.size = m_species.names.size();
    if (rank == 0)  std::cout << "    ";
    for(unsigned i =0; i < m_species.size; ++i){
      if (rank == 0)  std::cout << m_species.names[i] << " ";
      //m_species.indices[m_species.names[i]] = i;
      Trim(m_species.names[i]);
    }
    if (rank == 0)  std::cout << std::endl;
    
  }

  // get solid phases
  ///////////////////
  if (rank == 0)  std::cout << "Solid Phases" << std::endl;
  HierarchicalDataNode* phasesNode = hdn->GetChild("SolidPhases");
  {
    m_solidphases.names = phasesNode->GetStringVector("names", ",");
    m_solidphases.volumeFractions = phasesNode->GetAttributeVector<realT>("volume_fractions", ",");
    Trim(m_solidphases.names);
    m_solidphases.reactionRateFunctionNames = phasesNode->GetStringVector("reactionfunctions", ",");
    m_solidphases.size = m_solidphases.names.size();
    for(unsigned i =0; i < m_solidphases.size; ++i){
      if (rank == 0)  std::cout << m_solidphases.names[i] << " ";
    }

    Trim(m_solidphases.reactionRateFunctionNames);

    if (rank == 0)  std::cout << std::endl;
  }
  m_solidphases.reactionRateFunctionPtrs.resize(m_solidphases.size); 


  //GPChemistry::ReadXMLEquilibriumData(hdn,m_speciesIndexMap,m_equilibrium.equations,m_equilibrium.logK);
  
  
  // face set
  m_faceSetName = hdn->GetAttributeString("faceset");
  
  m_temperature = hdn->GetAttributeOrDefault("temperature",60);
  
  m_traceConcentration = hdn->GetAttributeOrDefault("trace_concentration",1e-32);

}



void GeochemicalReactionSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{

  // Unit scale
  /////////////
  {
    using namespace GPUnits;
    UnitManager& um = UnitManager::Instance(); 
    m_convertToMolPerL = um.ConvertTo("mol/l",1);
    m_convertFromMolPerL =  1.0/m_convertToMolPerL;
    m_convertFromMolPerLSMsqrd = um.Convert("mol/(s*l*m^2)"); // reaction rate per surface area conversion
  }


  m_faceSet = &(domain.m_feFaceManager.GetSet(m_faceSetName));
  
  ChemistryManager& chemistryManager = ChemistryManager::Instance();
  std::vector<GPChemistry::ChemicalReaction>& solidPhaseReactionList = chemistryManager.m_phases;

  m_reactionSpecies.insert(m_species.names.begin(),m_species.names.end());
  
  for(unsigned i =0; i < m_solidphases.size; ++i)
  { 
    std::vector<GPChemistry::ChemicalReaction>::iterator reactionItr = solidPhaseReactionList.begin();
    bool notFound = true;
    while( notFound && reactionItr != solidPhaseReactionList.end() )
    {
      notFound = !streq(m_solidphases.names[i],reactionItr->name);
      if(notFound) ++reactionItr;
    }
  	
    if(notFound)
    {
      std::cout << "Warning solid phase " << m_solidphases.names[i] 
                << " was not found in chemical database" << std::endl;
    }else{
      m_solidphases.reactions.push_back(*reactionItr);
      m_reactionSpecies.insert(reactionItr->species.begin(),reactionItr->species.end()); 
    }
  }

  BuildEquilibriumEquations();
  
  // get reaction functions
  FunctionManager& functionManager = FunctionManager::Instance();  
  
  for(unsigned i =0; i < m_solidphases.size; ++i){
  	m_solidphases.reactionRateFunctionPtrs[i] =
  	   &( functionManager.GetFunction( m_solidphases.reactionRateFunctionNames[i] ) );
  }
  
  // add element fields to faces
  for(unsigned i =0; i < m_elements.size; ++i){
    domain.m_feFaceManager.AddKeylessDataField<realT>(m_elements.names[i],true,true);
    domain.m_feFaceManager.AddKeylessDataField<realT>(m_elements.names[i]+"_ReactionRate", true,true);
  }
  // add solid phase reaction rates to faces
  for(unsigned i =0; i < m_solidphases.size; ++i){
    domain.m_feFaceManager.AddKeylessDataField<realT>(m_solidphases.names[i]+"_ReactionRate", true,true);
  }

  // equilibrium solver
  unsigned n = m_indexSpeciesMap.size();
  unsigned nk = m_equilibriumEquations.Dimension(0);
  unsigned na = n-nk;
  
  rArray2d& KK = m_equilibriumEquations;
  rArray1d pp(nk);
  for(unsigned i=0; i < nk; ++i){
    pp(i) = log(10)* m_LogK(i); // convert log_10 data to log_e
  }
   
  // element concentrations
  rArray2d AA(na,n);
  rArray1d bb(na,0.0); // eventually will be set to brine conc
  int ii = 0;
  for(unsigned i=0; i < na-1; ++i){
    std::string& el = m_elements.names[i];
    bool isInvalidEl = streq(el,"e") || streq(el,"H") || streq(el,"O");
    if(!isInvalidEl){
      for(unsigned j=0; j < n; ++j) AA(ii,j) = m_elements.concentrationEquations(i,j);
      ++ii;
    }
  }

  // charge balance
  for(unsigned j=0; j < n; ++j){
    AA(na-1,j) = m_valences(j);
  }
  bb(na-1) = 0.0;   


  // equilibrium equations
  m_brineSystem = ODRF::EquilibriumSystem(KK, pp, AA, bb,m_indexSpeciesMap,m_temperature);
    
}


void GeochemicalReactionSolver::BuildEquilibriumEquations()
{
  ChemistryManager& chemistryManager = ChemistryManager::Instance();

  // collect reactions
  chemistryManager.BuildReactionList(m_reactionSpecies, m_solutionReactions,m_solidPhaseReactions);
  

  // separate master species
  for(unsigned i =0; i < m_solutionReactions.size(); ++i)
  {
     m_solutionSpecies.insert(m_solutionReactions[i].species.begin(),
                              m_solutionReactions[i].species.end());
  }

  // remove H2O from list
  m_reactionSpecies.erase("H2O");
  m_solutionSpecies.erase("H2O");
  
  sArray1d speciesList;
  sArray1d masterSpeciesList;
  for(std::set<std::string>::iterator itr = m_solutionSpecies.begin(); 
      itr != m_solutionSpecies.end(); ++itr)
  {
    bool isMasterSpecies = false;
    for(unsigned i = 0; i < chemistryManager.m_solutionMasterSpecies.size() && !isMasterSpecies; ++i)
    {
    	isMasterSpecies = streq( chemistryManager.m_solutionMasterSpecies[i].species, *itr);
    }  	
      	
    if(isMasterSpecies )
    {
      m_solutionMasterSpecies.insert(*itr);
      masterSpeciesList.push_back(*itr);
    } else {
      speciesList.push_back(*itr);
    }
  }

  for(unsigned i =0; i < speciesList.size(); ++i)
  {
    m_speciesIndexMap[speciesList[i]] = i;
    m_indexSpeciesMap.push_back(speciesList[i]);
  }
  
  for(unsigned i =0; i < masterSpeciesList.size(); ++i)
  {
    m_speciesIndexMap[masterSpeciesList[i]] = i+speciesList.size();
    m_indexSpeciesMap.push_back(masterSpeciesList[i]);
  }

  // build reaction equations
  {
    unsigned m = m_solutionReactions.size();
    unsigned n = m_speciesIndexMap.size();
    rArray2d A(m,n); 
    rArray1d logK(m);
    for(unsigned i =0; i < m; ++i)
    {
      GPChemistry::ChemicalReaction& reaction = m_solutionReactions[i];
      for(unsigned jj = 0; jj < reaction.species.size(); ++jj)
      {
        if(!streq(reaction.species[jj],"H2O") ){
          unsigned j = m_speciesIndexMap[reaction.species[jj]];
          A(i,j) = reaction.stoich_coeffs[jj];
        }
      }
      logK(i) = reaction.EquilibriumConstant(m_temperature);
    }
  
    std::cout << std::endl;
    for(unsigned i =0; i < m; ++i)
    {
      for(unsigned j =0; j < n; ++j) std::cout << A(i,j) << " ";
      std::cout << " = " << logK(i) << std::endl;
    }

    // reduced row echelon form
    //std::vector<int> pivot_cols;
    //reducedRowEchelonForm(A, logK,pivot_cols);

    m_LogK = logK;
    m_equilibriumEquations = A;
  
    std::cout << std::endl;
    for(unsigned i =0; i < m; ++i)
    {
      for(unsigned j =0; j < n; ++j) std::cout << m_equilibriumEquations(i,j) << " ";
      std::cout << " = " << m_LogK(i) << std::endl;
    }
  }


  // build solid phase reactions
  {
    unsigned m = m_solidphases.reactions.size();
    unsigned n = m_speciesIndexMap.size();
    rArray2d A(m,n); 
    rArray1d logK(m);
    for(unsigned i =0; i < m; ++i)
    {
      GPChemistry::ChemicalReaction& reaction = m_solidphases.reactions[i];
      for(unsigned jj = 0; jj < reaction.species.size(); ++jj)
      {
        if(isMember(reaction.species[jj],m_speciesIndexMap) ){ // ignore solid phase & water
          unsigned j = m_speciesIndexMap[reaction.species[jj]];
          A(i,j) = reaction.stoich_coeffs[jj];
        }
      }
      logK(i) = reaction.EquilibriumConstant(m_temperature);
    }

    m_solidphases.logK = logK;
    m_solidphases.equilibriumEquations = A;
  }

    
  // elements
  unsigned totalNumSpecies = m_indexSpeciesMap.size();
  m_speciesElements.resize(totalNumSpecies);
  for(unsigned i=0; i < totalNumSpecies; ++i)
  {
    std::cout << m_indexSpeciesMap[i] <<std::endl; 
    GPChemistry::GetElementCountFromFormula(m_indexSpeciesMap[i], m_speciesElements[i]);
    
    std::map<std::string, realT>::iterator
      itr = m_speciesElements[i].begin(),
      iend = m_speciesElements[i].end();
    for(;itr != iend; ++itr){
      std::cout << "    " << itr->first << " " << itr->second << std::endl;
      if( !( streq(itr->first,"H") || streq(itr->first,"O") ) ) 
         m_elements.indexMap[itr->first]=0;
    }
  }
  std::cout<< std::endl << std::endl; 
 
  // record element indicies, names
  {  
    m_elements.size = m_elements.indexMap.size();
    m_elements.names.resize(m_elements.size);
    
    std::map<std::string, unsigned>::iterator 
      itr = m_elements.indexMap.begin(),
      iend = m_elements.indexMap.end();
    unsigned indx = 0;
    for(;itr!=iend;++itr,++indx){
      itr->second = indx;
      m_elements.names[indx] = itr->first;
    }
  }
  
  // build concentration equations
  m_elements.concentrationEquations = rArray2d(m_elements.size,totalNumSpecies);
  for(unsigned spIndx=0; spIndx < totalNumSpecies; ++spIndx)
  {
    std::map<std::string, realT>::iterator
      itr = m_speciesElements[spIndx].begin(),
      iend = m_speciesElements[spIndx].end();
    for(;itr != iend; ++itr)
    {
      const std::string& el = itr->first;
      realT count = itr->second;
      if( isMember(el,m_elements.indexMap) ){ // ignores H,O
        unsigned elIndx =  m_elements.indexMap[el];

        m_elements.concentrationEquations(elIndx,spIndx) = count;
      }
    }
  }

  // solid phase reaction elements
  {
    rArray2d A(m_solidphases.size,m_elements.size);
    for(unsigned i=0; i < m_solidphases.size; ++i)
    {
      std::cout << m_solidphases.names[i] <<std::endl; 
      std::string solidFormula = m_solidphases.reactions[i].species[0];
      std::cout << solidFormula <<std::endl; 
      std::map<std::string, realT> elMap;
      GPChemistry::GetElementCountFromFormula(solidFormula, elMap);
    
      std::map<std::string, realT>::iterator
        itr = elMap.begin(),
        iend = elMap.end();
      for(;itr != iend; ++itr){
        std::cout << "    " << itr->first << " " << itr->second << std::endl;
        if( !( streq(itr->first,"H") || streq(itr->first,"O") ) ) {
          unsigned index = m_elements.indexMap[itr->first];
          A(i,index) = itr->second;
        }
      }
    }
    std::cout<< std::endl << std::endl;

    m_solidphases.reactionElements=A;

  }
  
  // get species valences
  m_valences.resize(totalNumSpecies);
  for(unsigned spIndx=0; spIndx < totalNumSpecies; ++spIndx){
  	m_valences[spIndx] = GPChemistry::GetValence(m_indexSpeciesMap[spIndx]);
  }
  
  
  for(unsigned i=0; i < m_elements.size; ++i){
    std::cout << m_elements.names[i] << ":\t ";
    for(unsigned j=0; j < totalNumSpecies; ++j){
      std::cout << m_elements.concentrationEquations(i,j) << " " ;
    }
    std::cout << std::endl;
  } 	
  std::cout << std::endl;
  
  
  std::cout << "Solution species" << std::endl;
  for(unsigned i=0; i < m_solutionReactions.size(); ++i)
  {
    std::cout << m_solutionReactions[i].name <<": "<< m_solutionReactions[i].formula <<std::endl; 
  }
  std::cout<< std::endl << std::endl; 


  std::cout << "Solid species" << std::endl;
  for(unsigned i=0; i < m_solidphases.size; ++i)
  {
    std::cout << m_solidphases.names[i] <<": "<< m_solidphases.reactions[i].formula <<std::endl; 
  }
  std::cout<< std::endl << std::endl; 
  
  std::cout << "Valences" << std::endl;
  for(unsigned spIndx=0; spIndx < totalNumSpecies; ++spIndx){
  	std::cout << "    " <<  m_indexSpeciesMap[spIndx] << " "<< m_valences[spIndx] << "\n";
  }
  std::cout<< std::endl << std::endl; 
  
  
}


void GeochemicalReactionSolver::TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber, PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator ){

  FunctionManager& functionManager = FunctionManager::Instance(); 

  unsigned numSpecies =  m_indexSpeciesMap.size(); 
  
  rArray1d solidphaseReactionRates(m_solidphases.size,0.0);
  rArray1d elementReactionRates(m_elements.size,0.0);  
  rArray1d activities(numSpecies,0.0);
  rArray1d molalities(numSpecies*2,m_traceConcentration); // initial guess
  for(unsigned i =0; i < numSpecies; ++i) molalities(i+numSpecies) = log(std::max(molalities(i),1e-64));

  std::vector< rArray1d* > elementFieldPtrs(m_elements.size );
  std::vector< rArray1d* > elementReactionRateFieldPtrs(m_elements.size );
  for(unsigned i =0; i < m_elements.size; ++i){
    elementFieldPtrs[i] = &(domain.m_feFaceManager.GetFieldData<realT>( m_elements.names[i] ));
    elementReactionRateFieldPtrs[i] 
           = &(domain.m_feFaceManager.GetFieldData<realT>( m_elements.names[i] +"_ReactionRate" ));
  }

  std::vector< rArray1d* > solidReactionRateFieldPtrs(m_solidphases.size );
  for(unsigned i =0; i < m_solidphases.size; ++i){
    solidReactionRateFieldPtrs[i] 
           = &(domain.m_feFaceManager.GetFieldData<realT>( m_solidphases.names[i] +"_ReactionRate" ));
  }
  

  // reaction rate functions are of the form
  // f(SurfaceArea, Ion activity product (Q) , list of activities );
  //
  // we need reaction rate per unit surface area, so the
  // surface area supplied to the function is replaced by
  // the mineral volume fraction. 
  //
  unsigned nVars = m_solidphases.reactionRateFunctionVariables.size();
  std::vector<realT> x(nVars+2); // function buffer
  std::vector<unsigned> packIndxs(nVars);
  for(localIndex i = 0; i < nVars; ++i){
    std::string& var = m_solidphases.reactionRateFunctionVariables[i];
    packIndxs[i] = m_speciesIndexMap[var];
  }
    
  // loop over set
  lSet::const_iterator 
    fi=m_faceSet->begin(),
    fiend=m_faceSet->end();
  for(;fi != fiend; ++fi){
  	
    const localIndex& fc = *fi;
    
    // update element concentrations
    rArray1d& b = m_brineSystem.Get_b();
    for(unsigned i =0; i < m_elements.size; ++i){
      b[i] = std::max( (*elementFieldPtrs[i])[fc], m_traceConcentration );
      b[i] *= m_convertToMolPerL;
    }

    // get equilibrium
    const unsigned maxNumItrs = 300;
    const realT tolx = std::numeric_limits<realT>::epsilon();

    NewtonRaphson( m_brineSystem, molalities, maxNumItrs, tolx);
    
    // calculate activities
    activities = molalities; // fixme

    // calculate reaction rates
    {
	
      // pack function activities into x
      for(localIndex i = 0; i < nVars; ++i){
      	/* std::string& var = m_solidphases.variables[i];
      	   unsigned index = m_speciesIndexMap[var]; */
  	// pack field values into x
        x[i+2] = activities[packIndxs[i]];
      }
    
      for(unsigned i =0; i <m_solidphases.size ; ++i){

        // effective surface area = volume fraction
        x[0] = m_solidphases.volumeFractions[i];
      
        // ion activity product
  	double Q = 1.0; 
  	for(unsigned j =0; j < numSpecies ; ++j){
  	   realT stoich = m_solidphases.equilibriumEquations(i,j);
  	   if(stoich != 0.0) Q *= pow(activities[j],stoich);
  	}
      
        x[1] = Q;
      	
        Function& func = *(m_solidphases.reactionRateFunctionPtrs[i]);	
        solidphaseReactionRates[i] = func(x[0])*m_convertFromMolPerLSMsqrd;
      }
    }
    
    // apply reaction rates to net dissolved element fields
    for(unsigned i =0; i< m_solidphases.size; ++i){
      for(unsigned j =0; j< m_elements.size; ++j){
        realT coeff = m_solidphases.reactionElements(i,j);
        if(coeff > 0.0)
          elementReactionRates[j] += coeff*solidphaseReactionRates[i];
      }
    }

    // copy reaction rates to field
    for(unsigned i =0; i< m_elements.size; ++i){
       (*elementReactionRateFieldPtrs[i])[fc] = elementReactionRates[i];
    }

    // copy solid phase reaction rates to field
    for(unsigned i =0; i< m_solidphases.size; ++i){
       (*solidReactionRateFieldPtrs[i])[fc] = solidphaseReactionRates[i];
    }
    
    
  }

};


REGISTER_SOLVER( GeochemicalReactionSolver )

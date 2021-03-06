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
//  LLNL-CODE-656616
//  GEOS-CORE, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GEOS-CORE. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//
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
 * @file UnitManager.h
 * @author walsh24
 * @date July 27, 2011
 */

#ifndef UNITMANAGER_H_
#define UNITMANAGER_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "Utilities/Functions.h"
#include <string>
#include <cmath>

namespace TICPP { class HierarchicalDataNode; }

/// Units in the GPUnits namespace are expressed in m,g,s,K,mol
/// For unit conversion use the UnitManager class. 
namespace GPUnits
{
  enum basic_dimensions {
     UNIT_LENGTH, UNIT_MASS, UNIT_TIME, UNIT_TEMPERATURE, UNIT_MOLE, NUM_UNITS
  };
  
  // SI prefixes
  const realT giga = 1e9;
  const realT mega = 1e6;
  const realT kilo = 1e3; 
  const realT deca = 10.0; 
  const realT deci = 0.1; 
  const realT centi = 1e-2; 
  const realT milli = 1e-3; 
  const realT micro = 1e-6; 
  const realT nano = 1e-9; 

  // length units
  const realT m = 1.0;
  const realT cm = centi*m;
  const realT mm = milli*m;
  const realT km = kilo*m;
  const realT um = micro*m; // micrometers
  const realT nm = nano*m;
  const realT in = 2.54e-2*m;
  const realT ft = 0.3048*m; //12.0*in
  const realT yd = 36.0*in;
  const realT mi = 1760.0*yd;

  // area/permeability units
  const realT darcy = 9.869233e-13*m*m;
  
  
  // mass units
  const realT g = 1.0;
  const realT mg = milli*g;
  const realT kg = kilo*g;
  const realT lb = 453.59237*g;
  const realT slug = 14.593903*kg; // unit of mass used in lbf calculations
  
  // time units
  const realT s = 1.0; 
  const realT ms = milli*s;
  const realT us = micro*s; // microseconds
  const realT min = 60*s; 
  const realT hour = 60*min; 
  const realT day = 24*hour; 
  const realT year = 365.25*day;// Julian year
  const realT shake = 1e-8*s; // shake = 10 nano seconds;
  
  // temperature units
  const realT K = 1.0;
  
  // amount units
  const realT mol = 1.0;
  const realT mmol = milli*mol;

  struct UnitData{
	  realT value;
	  std::string name;
  };
 
}

/**
 * @author walsh24
 * 
 * @brief Singleton class to control unit conversion and set default units. 
 * 
 **/
class UnitManager{
	
  public:
    static UnitManager& Instance()
    {
      static UnitManager theUnitManager;
      return theUnitManager;
    };
    void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;
    
    /// Convert string expression to base units
    realT Convert(const std::string& quantityString);
    /// Convert quantity and string units to base units
    realT Convert(realT value, std::string unitString);
    /// Convert quantity with explicitly defined units to base units (faster than string methods)
    realT Convert(realT value, realT lengthUnit, int lengthdim, 
                               realT massUnit= 1.0, int massdim = 0,
                               realT timeUnit= 1.0, int timedim = 0,
                               realT tempUnit= 1.0, int tempdim = 0,
                               realT moleUnit= 1.0, int moledim = 0);
    
    
    /// Convert quantity from base units to units given in a string
    realT ConvertTo(const std::string& unitString, realT value);
    /// Convert quantity from base units to explicitly defined units  (faster than string methods)
    realT ConvertTo(realT lengthUnit, int lengthdim, 
                    realT massUnit, int massdim,
                    realT timeUnit, int timedim,
                    realT tempUnit, int tempdim,
                    realT moleUnit, int moledim, realT value);
    // overloaded functions to allow ConvertTo to mirror the form of Convert
    realT ConvertTo(realT lengthUnit, int lengthdim, realT massUnit, int massdim, realT timeUnit, int timedim,
                    realT tempUnit, int tempdim, realT value)
    { 
      return ConvertTo(lengthUnit, lengthdim, massUnit, massdim, timeUnit, timedim,tempUnit, tempdim, 1.0, 0, value);
    };
    realT ConvertTo(realT lengthUnit, int lengthdim, realT massUnit, int massdim, realT timeUnit, int timedim, realT value)
    { 
      return ConvertTo(lengthUnit, lengthdim, massUnit, massdim, timeUnit, timedim,1.0, 0, 1.0, 0, value);
    };
    realT ConvertTo(realT lengthUnit, int lengthdim, realT massUnit, int massdim, realT value)
    { 
      return ConvertTo(lengthUnit, lengthdim, massUnit, massdim, 1.0, 0,1.0, 0, 1.0, 0, value);
    };
    realT ConvertTo(realT lengthUnit, int lengthdim, realT value)
    { 
      return ConvertTo(lengthUnit, lengthdim, 1.0, 0, 1.0, 0,1.0, 0, 1.0, 0, value);
    };


    void ReportAllUnits(void){
      std::map<std::string,GPUnits::UnitData>::iterator itr =  m_unitValueMap.begin();
	  std::map<std::string,GPUnits::UnitData>::iterator iend = m_unitValueMap.end() ;
      std::cout << std::endl;
      std::cout << "Symbol" << "\t \t" << "Description"<< "\n";
      std::cout << "------" << "\t \t" << "-----------"<< "\n";
      for(;itr != iend; ++itr){
        std::cout << itr->first << "\t \t" << itr->second.name << "\n";
      }
      std::cout << std::endl;
    };
 
  private:


 
    /**
     *  @brief Constructor
     * 
     *  NB. The units are set to SI units by default, 
     *  Within the class, the base units are expressed in m,g,s,K,mol.
     * 
     **/
    UnitManager():
      m_fParser()
    {
      using namespace GPUnits;
    	// SI units - m kg s K mol
      m_baseUnits[UNIT_LENGTH]= m;
      m_baseUnits[UNIT_MASS]= kg;
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
      UpdateParser();
    };

    ~UnitManager() {};
    UnitManager( const UnitManager& );
    UnitManager& operator=( const UnitManager& );

    realT Units( realT length= 0.0,
                 realT mass= 0.0,
                 realT time= 0.0,
                 realT temp  = 0.0,
                 realT mole =0.0)
    {
       realT rv = 1.0;
       if( !isZero(length) ) rv /= pow(m_baseUnits[0],length);
       if( !isZero(mass) )   rv /= pow(m_baseUnits[1],mass);
       if( !isZero(time) )   rv /= pow(m_baseUnits[2],time);
       if( !isZero(mole) )   rv /= pow(m_baseUnits[3],mole);
       return rv;
    };
    void UpdateParser(void);
    
    /// Add a new unit symbol to the unit manager
    void DefineUnit(const std::string& name,const std::string& value,const std::string& description="");
    void DefineUnit(const std::string& name,const realT& value,const std::string& description="");


    realT m_baseUnits[GPUnits::NUM_UNITS]; // problem default units given in m,g,s,K,mol. 
    FunctionParser m_fParser;
    std::map<std::string,GPUnits::UnitData> m_unitValueMap;
};

//////////////////////////////////////////////////////////////////////////////////////////

/**
 * @author walsh24
 * @brief Converts a quantity in the default units to another set of units.
 * 
 * eg.
 *  UnitManager& theUnitManager = UnitManager::Instance();
 *  realT gbase = theUnitManager.Convert("9.8 m/s^2");  // in base units
 *  realT g = theUnitManager.ConvertTo("km/s^2",gbase);  // g = 0.0098
 * 
 **/
 inline
realT UnitManager::ConvertTo(const std::string& unitString,realT value){
	realT invVal = Convert(1.0, unitString);
	return value/invVal;
}

/**
 * @author walsh24
 * @brief Converts a quantity in the default units to another set of units.
 * 
 * eg.
 *  UnitManager& theUnitManager = UnitManager::Instance();
 *  realT gbase = theUnitManager.Convert("9.8 m/s^2");  // in base units
 *  using namespace GPUnits;
 *  realT g = theUnitManager.ConvertTo(km,1,s,-2,gbase);  // g = 0.0098
 * 
 **/
 inline
realT UnitManager::ConvertTo(realT lengthUnit, int lengthdim, 
                             realT massUnit, int massdim, realT timeUnit, int timedim,
                             realT tempUnit, int tempdim, realT moleUnit, int moledim, realT value){    
   
   realT invVal = Convert(1.0, lengthUnit, lengthdim, massUnit, massdim, timeUnit, timedim,
                              tempUnit, tempdim, moleUnit, moledim);
   return value/invVal;                     	
}





#endif /* UNITMANAGER_H_ */

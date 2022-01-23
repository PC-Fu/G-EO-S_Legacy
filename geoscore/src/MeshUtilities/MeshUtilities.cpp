/*
 * MeshUtilities.cpp
 *
 *  Created on: Dec 5, 2012
 *      Author: settgast1
 */

#include "MeshUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "SimpleGeometricObjects.h"

MeshUtilities::MeshUtilities()
{
  // TODO Auto-generated constructor stub

}

MeshUtilities::~MeshUtilities()
{
  // TODO Auto-generated destructor stub
}



void MeshUtilities::GenerateNodesets( TICPP::HierarchicalDataNode& hdn,
                                      NodeManagerT& nodeManager )
{

  std::map< std::string, lSet >& sets = nodeManager.m_Sets;
  Array1dT<R1Tensor>& X = *(nodeManager.m_refposition);

  for(TICPP::HierarchicalDataNode* hdnNode = hdn.Next(true); hdnNode; hdnNode = hdn.Next() )
  {
    std::string header = "Nodeset";
    if( !( header.compare(hdnNode->Heading() ) ) )
    {
      /*
      SimpleGeometricObjectBase::Types type = SimpleGeometricObjectBase::IntToType(hdnNode->GetAttributeValue<int>("type") );
      */

      std::string name = hdnNode->GetAttributeString("name");
      lSet& set = sets[name];

      SimpleGeometricObjectBase* object;

      // new allocation method
      TICPP::HierarchicalDataNode* geometryNode = hdnNode->GetChild("Geometry");
      if(geometryNode){
    	// treating geometry node as a union boolean geometry allows multiple objects to be defined within it:
    	  /**
    	    <Geometry>
    	        <Box  ... />
    	        <Sphere .. />
    	        <Intersection>
    	          <Cylinder .. />
    	          <Not>
    	            <Cylinder .. />
    	          </Not>
    	        </Intersection>
    	    <Geometry>
    	    **/
        object = SimpleGeometricObjectBase::Allocate(SimpleGeometricObjectBase::unionGeometry);
        object->ReadXML( *geometryNode );
      } else {
        // old allocation method for backwards compatability
        std::string geometricObjectTypeStr = hdnNode->GetAttributeStringOrDefault("type", "Box");
        SimpleGeometricObjectBase::Types type = fromString<SimpleGeometricObjectBase::Types>(geometricObjectTypeStr);
        object = SimpleGeometricObjectBase::Allocate( type );

        object->ReadXML( *hdnNode );
      }

      for( localIndex a=0 ; a<X.size() ; ++a )
      {
        if( object->IsCoordInObject( X[a] ) )
        {
          set.insert(a);
        }
      }


      SimpleGeometricObjectBase::Deallocate( object );
    }
  }
}

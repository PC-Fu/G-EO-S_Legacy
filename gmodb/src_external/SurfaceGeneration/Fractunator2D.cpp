/**
 * @file Fractunator2D.cpp, modified based on FractunatorBase.cpp
 * @author fu4
 * @date Oct. 29, 2012
 */

#include "Fractunator2D.h"
//#include "Common/typedefs.h"
#include <limits.h>
#include <mpi.h>
#include <math.h>
//#include "IO/BinStream.h"


#include "MPI_Communications/SpatialPartition.h"
#include "DataStructures/VectorFields/ElementRegionT.h"
#include "Utilities/GeometryUtilities.h"

#include "SurfaceGeneration/FractunatorFactory.h"







Fractunator2D::Fractunator2D():
FractunatorBase(),
m_maxKinkAngle(140.0),
m_kinkStrength(0.0)

{
  // TODO Auto-generated constructor stub

}

Fractunator2D::~Fractunator2D()
{
  // TODO Auto-generated destructor stub
}


void Fractunator2D::RegisterFieldsAndMaps( NodeManagerT& nodeManager,
                                           EdgeManagerT& edgeManager,
                                           FaceManagerT& faceManager )
{

  // 2D fractunator is not calling the base function.
  // At some point we will want to clean the relationship between 2D, 3, and base and call the base function here too.

  // ruptureState is equivalent to flagFrEdge in FracFlow
  faceManager.AddKeylessDataField<int>( "ruptureState",true, true );

  nodeManager.AddKeylessDataField<int>("numberOfRupturedFaces",true,true);

  nodeManager.AddKeylessDataField<int>("isDetachedFromSolidMesh", true, false );
  iArray1d& isDetachedFromSolidMesh = nodeManager.GetFieldData<int>("isDetachedFromSolidMesh");
  isDetachedFromSolidMesh = 0;

  nodeManager.AddKeylessDataField<int>("nodeSplitability", true, true);
  // = -1 if we haven't evaluated this node yet.
  // = 0, this node is known to be unsplitable (either an interior node or the kink angle is too large.
  // = 1, this node is a fracture front tip.  We might need to check the saturation of the cell connected to it.
  // = 2, this node is at a kink with a kink angle that is small enough.
  // We only evaluate this variable when we split a node.
  //  iArray1d& nodeSplitability = nodeManager.GetFieldData<int>("nodeSplitability");
  isDetachedFromSolidMesh = 0;

  nodeManager.AddKeylessDataField<realT>("kinkAngle", true, false );
  nodeManager.AddKeylessDataField<realT>("SIF_I", false, true );
  nodeManager.AddKeylessDataField<realT>("SIF_II", false, true );
  faceManager.AddKeylessDataField<realT>("SIFonFace", true, true );


  nodeManager.AddKeylessDataField<R1Tensor>("displacement0", true, false );
  nodeManager.AddKeylessDataField<R1Tensor>("netDisplacement", false, true );


  //nodeManager.AddMap< UnorderedVariableOneToManyRelation >("usedFaces");
  nodeManager.AddMap<OrderedVariableOneToManyRelation>("childIndices");
  nodeManager.AddMap<OneToOneRelation>("parentIndex");
  OneToOneRelation& parentIndexNodes = nodeManager.GetOneToOneMap( "parentIndex" );
  parentIndexNodes = LOCALINDEX_MAX;

  //  I suspect that we don't need to track child/parent of edges.  We will see.
  //  edgeManager.AddMap<OrderedVariableOneToManyRelation>("childIndices");
  //  edgeManager.AddMap<OneToOneRelation>("parentIndex");
  //  OneToOneRelation& parentIndexEdge = edgeManager.GetOneToOneMap( "parentIndex" );
  //  parentIndexEdge = LOCALINDEX_MAX;

  faceManager.AddMap<OrderedVariableOneToManyRelation>("childIndices");
  faceManager.AddMap<OneToOneRelation>("parentIndex");
  OneToOneRelation& parentIndexFace = faceManager.GetOneToOneMap( "parentIndex" );
  parentIndexFace = LOCALINDEX_MAX;

  faceManager.AddKeylessDataField<realT>( "separationCoeff",true, true );
  faceManager.AddKeylessDataField<realT>( "stressNOnFace",true, true );
  faceManager.AddKeylessDataField<R1Tensor>( "stressTOnFace",true, true );
  rArray1d& stressNOnFace = faceManager.GetFieldData<realT>("stressNOnFace");
  stressNOnFace = -std::numeric_limits<realT>::max();



  // TODO Is there a better way handling this?  The field should be registered by ParallelPlateFlowSolverBase.  We have to register here is we do fracturing without a flow solver.
  edgeManager.AddMap< UnorderedVariableOneToManyRelation >( "edgeToFlowFaces");


  faceManager.AddKeylessDataField<realT>( "birthTime",true, true );
  faceManager.AddKeylessDataField<int>("isSeparable", true, true );


  rArray1d* toughnessPointer = faceManager.GetFieldDataPointer<realT>("faceToughness");
  if (toughnessPointer == NULL)
  {
    faceManager.AddKeylessDataField<realT>("faceToughness", true, true);
    m_tounessSetByInitialCondition = false;
  }
  else
  {
    m_tounessSetByInitialCondition = true;
  }

  rArray1d* failStressPointer = faceManager.GetFieldDataPointer<realT>("faceFailStress");
  if (failStressPointer == NULL)
  {
    if (m_failCriterion == 0 || m_failCriterion == 2)
    {
      faceManager.AddKeylessDataField<realT>("faceFailStress", true, true);
    }
    m_failStressSetByInitialCondition = false;
  }
  else
  {
    m_failStressSetByInitialCondition = true;
  }

}


void Fractunator2D::Initialize( NodeManagerT& nodeManager,
                                EdgeManagerT& edgeManager,
                                FaceManagerT& faceManager,
                                ElementManagerT& elementManager)
{

  FractunatorBase::Initialize( nodeManager, edgeManager, faceManager, elementManager );

  //Move reference location to create preexisting stress field.
  Array1dT<R1Tensor>& u = nodeManager.GetFieldData<FieldInfo::displacement>();
  Array1dT<R1Tensor>& x0 = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  Array1dT<R1Tensor>& du = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();

  Array1dT<R1Tensor>& u0 = nodeManager.GetFieldData<R1Tensor>("displacement0");


  //Move the ref position of nodes to generate the desired stress field.
  // For now, we assume plane strain
  // This only works for uniform material properties
  {
    const ElementRegionT& elemRegion = (elementManager.m_ElementRegions.begin())->second;
    const MaterialBaseParameterData& parameter = *(elemRegion.m_mat->ParameterData(0));
    const realT E = parameter.E;
    const realT v = parameter.Nu;
    R1Tensor s = m_insituStress_2D;
    s *= 1+v;
    s /= E;

    const realT strainxx = (1-v) * s[0] - v * s[1] ;
    const realT strainyy = -v * s[0] + (1-v) * s[1];
    const realT strainxy = s[2];
    const realT strainyx = strainxy;

    for (localIndex i = 0; i != nodeManager.DataLengths(); ++i)
    {
      R1Tensor x_current = x0[i];
      x0[i][0] = x_current[0] / (1 + strainxx);
      x0[i][1] = x_current[1] / (1 + strainyy);
      x0[i][0] = (x0[i][0] - x0[i][1] * strainxy) / (1 - strainxy * strainyx);
      x0[i][1] = (x0[i][1] - x0[i][0] * strainyx) / (1 - strainxy * strainyx);

      u[i] = x_current;
      u[i] -= x0[i];
      u0[i] = u[i];
      du[i] = u[i];
    }
  }


  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegionIter = elementManager.m_ElementRegions.begin() ;
      elementRegionIter != elementManager.m_ElementRegions.end() ;
      ++elementRegionIter )
  {
    ElementRegionT& elementRegion = elementRegionIter->second;

    elementRegion.CalculateVelocityGradients(nodeManager);

    elementRegion.MaterialUpdate(0.0);

  }




  //  std::map< std::string, lSet >::const_iterator setMap = nodeManager.m_Sets.find( m_separableSet );
  //  iArray1d& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
  //  isSeparable = 1;
  //
  //
  //  // process nodes on the interior
  //  if( setMap!=nodeManager.m_Sets.end() )
  //  {
  //    isSeparable = 0;
  //    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
  //    {
  //      if( setMap->second.count(a)==1 )
  //        isSeparable[a] = 1;
  //    }
  //  }
  //
  //  rArray1d& maxCohesiveTraction = faceManager.GetFieldData<realT>("maxCohesiveTraction");
  //  maxCohesiveTraction = m_failstress;
  //
  //
  //  rArray1d& separationCoeff = faceManager.GetFieldData<realT>("separationCoeff");
  //  separationCoeff = 0.0;

}


void Fractunator2D::PreexistingFracture2D( NodeManagerT& nodeManager,
                                           EdgeManagerT& edgeManager,
                                           FaceManagerT& faceManager,
                                           ExternalFaceManagerT& externalFaceManager,
                                           ElementManagerT& elementManager,
                                           SpatialPartition& partition,
                                           const bool prefrac)
{


  Array1dT<lSet> nodesToRupturedFaces;
  Array1dT<lSet> edgesToRupturedFaces;

  // We call this just to calculate stress on faces. All rupture state should be zero after this.
  UpdateRuptureStates( nodeManager,
                       edgeManager,
                       faceManager,
                       elementManager,
                       nodesToRupturedFaces,
                       edgesToRupturedFaces,
                       prefrac );

  //MarkPreexistingFractures(faceManager, nodeManager, partition);

  CreatePreexistingFractures(nodeManager, edgeManager, faceManager, externalFaceManager, elementManager, partition, prefrac);


  rArray1d& kinkAngle = nodeManager.GetFieldData<realT>("kinkAngle");


  for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )

  {
    kinkAngle[a] = CalculateKinkAngle(a, nodeManager, faceManager);
  }


}



void Fractunator2D::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  FractunatorBase::ReadXML(hdn);

  m_failCriterion = hdn.GetAttributeOrDefault<int>("failCriterion",1);
  m_maxKinkAngle = hdn.GetAttributeOrDefault<realT>("maxKinkAngle",0.0);
  m_kinkStrength = hdn.GetAttributeOrDefault<realT>("kinkStrength",2.0e100);

  std::string temp = hdn.GetAttributeString("insitu_Stress");
  if( !temp.empty() )
  {
    m_insituStress_2D.StrVal( temp );
  }
  else
  {
    m_insituStress_2D = 0.0;
  }


  rArray1d x1_PreFrac, y1_PreFrac, z1_PreFrac, x2_PreFrac, y2_PreFrac, z2_PreFrac;
  x1_PreFrac = hdn.GetAttributeVector<realT>("x1_PreFrac", ",");
  y1_PreFrac = hdn.GetAttributeVector<realT>("y1_PreFrac", ",");
  z1_PreFrac = hdn.GetAttributeVector<realT>("z1_PreFrac", ",");
  x2_PreFrac = hdn.GetAttributeVector<realT>("x2_PreFrac", ",");
  y2_PreFrac = hdn.GetAttributeVector<realT>("y2_PreFrac", ",");
  z2_PreFrac = hdn.GetAttributeVector<realT>("z2_PreFrac", ",");

  m_x0_PreFrac.resize(x1_PreFrac.size());
  m_x1_PreFrac.resize(x2_PreFrac.size());

  for (localIndex iPF = 0; iPF < x1_PreFrac.size() ; ++iPF)
  {
    m_x0_PreFrac[iPF][0] = x1_PreFrac[iPF];
    m_x0_PreFrac[iPF][1] = y1_PreFrac[iPF];
    m_x0_PreFrac[iPF][2] = z1_PreFrac[iPF];
    m_x1_PreFrac[iPF][0] = x2_PreFrac[iPF];
    m_x1_PreFrac[iPF][1] = y2_PreFrac[iPF];
    m_x1_PreFrac[iPF][2] = z2_PreFrac[iPF];
  }


}




int Fractunator2D::SeparationDriver( NodeManagerT& nodeManager,
                                     EdgeManagerT& edgeManager,
                                     FaceManagerT& faceManager,
                                     ExternalFaceManagerT& externalFaceManager,
                                     ElementManagerT& elementManager,
                                     SpatialPartition& partition,
                                     const bool prefrac,
                                     const realT time)
{
  int rval = 0;
  Array1dT<R1Tensor>& displacement = nodeManager.GetFieldData<FieldInfo::displacement> ();
  const Array1dT<int>& isDetachedFromSolidMesh = nodeManager.GetFieldData<int> ("isDetachedFromSolidMesh");

  for (localIndex i = 0; i != nodeManager.DataLengths(); ++i)
  {
    if (isDetachedFromSolidMesh[i] ==1 && nodeManager.m_childIndices[i].size() == 2)
    {
      displacement[i] = displacement[nodeManager.m_childIndices[i][0]];
      displacement[i] += displacement[nodeManager.m_childIndices[i][1]];
      displacement[i] /=2;
    }
  }


  if (!prefrac)
  {
    if (m_failCriterion >0 )  // Stress intensity factor based criterion and mixed criterion.
    {
      if (m_failCriterion == 1)
      {
        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, std::numeric_limits<realT>::max()); //this->m_failstress );
      }
      else
      {
        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, m_failstress);
      }

      rArray1d& SIFonFace = faceManager.GetFieldData<realT>("SIFonFace");
      SIFonFace = std::numeric_limits<double>::min();

      if (m_failCriterion == 1)
      {
        iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
        for (localIndex a=0; a<faceManager.DataLengths(); ++a)
        {
          if (faceManager.m_toElementsRelation[a].size() == 2) ruptureState[a]=0;
        }
      }

      IdentifyRupturedFaces( nodeManager,
                             edgeManager,
                             faceManager,
                             elementManager,
                             partition,
                             prefrac );
    }
    else
    {
      // We cannot call the fractunator's UpdateRuptureStates function as in 3D.
      // Some extra stuff in that function does not apply to 2D.
      faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, m_failstress);
    }
  }
  else  // In the prefrac call, we need this to get the stressNOnFace, which will be used in the initialization of contacts for preexisting fractures.
  {
    for (localIndex kf = 0; kf < faceManager.DataLengths(); ++kf) faceManager.CalculateStressOnFace(elementManager, nodeManager, kf);
  }



  int rank ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  const iArray1d& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
  const iArray1d& layersFromDomainBoundary = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");



  // process nodes on the interior
  {
    ModifiedObjectLists modifiedObjects;
    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
    {
      //      const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);
      if( layersFromDomainBoundary[a]>1 &&
          (isSeparable[a] || prefrac)&&
          isNodeGhost[a]<0 &&
          nodeManager.m_toElementsRelation[a].size()>1 &&
          isDetachedFromSolidMesh[a] == 0 )
      {
        rval += EvaluateAndSplitNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, elementManager, modifiedObjects ) ;
      }
    }
    MarkBirthTime(faceManager, modifiedObjects, time);
  }


  for( int color=0 ; color<partition.NumColor() ; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( partition.Color() == color )
    {

      // process "near-boundary" nodes
      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
      {

        //        const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);

        if( layersFromDomainBoundary[a]<=1 &&
            (isSeparable[a] || prefrac) &&
            isNodeGhost[a]<0 &&
            nodeManager.m_toElementsRelation[a].size()>1 &&
            isDetachedFromSolidMesh[a] == 0 )
 //           nodesToRupturedFaces[a].size()>0 )
        {
          rval += EvaluateAndSplitNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, elementManager, modifiedObjects ) ;
        }
      }

      MarkBirthTime(faceManager, modifiedObjects, time);
    }
    partition.ModifyGhostsAndNeighborLists( modifiedObjects );
  } // color loop




//
//
//
//
//  const iArray1d& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();
//  rArray1d& kinkAngle = nodeManager.GetFieldData<realT>("kinkAngle");
//  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
//
//
//
//
//  const iArray1d& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
//  const iArray1d& layersFromDomainBoundary = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
//
//  Array1dT<lSet> nodesToRupturedFaces;
//  Array1dT<lSet> edgesToRupturedFaces;
//
//
//  UpdateRuptureStates( nodeManager,
//                       edgeManager,
//                       faceManager,
//                       elementManager,
//                       nodesToRupturedFaces,
//                       edgesToRupturedFaces,
//                       prefrac );
//
//  if (m_failCriterion == 1)
//  {
//    ruptureState = 0.0;
//  }
//
//
//  // process nodes on the interior
//  const localIndex  nNodes = nodeManager.DataLengths();  // By using nNodes instead of the length of nodeManager, we don't split child nodes that are created in the current step
//  {
//    ModifiedObjectLists modifiedObjects;
//    for( localIndex a=0 ; a<nNodes ; ++a )
//    {
//      if( layersFromDomainBoundary[a]>1 &&
//          isSeparable[a] &&
//          isNodeGhost[a]<0 &&
//          nodeManager.m_toElementsRelation[a].size() >= 2 &&
//          isDetachedFromSolidMesh[a] == 0)
//      {
//        if (CheckAndSplitNode( a,
//                           nodeManager,
//                           edgeManager,
//                           faceManager,
//                           externalFaceManager,
//                           elementManager,
//                           modifiedObjects,
//                           prefrac))
//        {
//          rval++;
//          std::cout << "rank " << rank << " split tip " << a << std::endl;
//        }
//
//      }
//    }
//
//    MarkBirthTime(faceManager, modifiedObjects, time);
//  }
//
//  // process "near-boundary" nodes
//  for( int color=0 ; color<partition.NumColor() ; ++color )
//  {
//    ModifiedObjectLists modifiedObjects;
//    if( partition.Color() == color )
//    {
//      for( localIndex a=0 ; a<nNodes ; ++a )
//      {
//        if( layersFromDomainBoundary[a]<=1 &&
//            isSeparable[a] &&
//            isNodeGhost[a]<0 &&
//            nodeManager.m_toElementsRelation[a].size() >= 2 &&
//            isDetachedFromSolidMesh[a] == 0)
//        {
//          if (CheckAndSplitNode( a,
//                                 nodeManager,
//                                 edgeManager,
//                                 faceManager,
//                                 externalFaceManager,
//                                 elementManager,
//                                 modifiedObjects,
//                                 prefrac) )
//          {
//            rval++;
//            std::cout << "rank " << rank << " split tip " << a << std::endl;
//          }
//        }
//      } // Node loop
//
//      MarkBirthTime(faceManager, modifiedObjects, time);
//
//
//      partition.ModifyGhostsAndNeighborLists( modifiedObjects );
//
//    }  // Color loop
//  }


  // We do this for all nodes for now.  There might be a smarter way based on the list of changed nodes passed from neighbors.
  rArray1d& kinkAngle = nodeManager.GetFieldData<realT>("kinkAngle");

  for( localIndex iNd = 0; iNd < nodeManager.DataLengths(); ++iNd )
  {
    if (isNodeGhost[iNd] < 0)
    {
      kinkAngle[iNd] = CalculateKinkAngle(iNd, nodeManager, faceManager);
    }
  }

  return rval;
}

void Fractunator2D::IdentifyRupturedFaces( NodeManagerT& nodeManager,
                                          EdgeManagerT& edgeManager,
                                          FaceManagerT& faceManager,
                                          ElementManagerT& elementManager,
                                          SpatialPartition& partition,
                                          const bool prefrac  )
{
  const iArray1d& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
  const iArray1d& layersFromDomainBoundary = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  const Array1dT<int>& isDetachedFromSolidMesh = nodeManager.GetFieldData<int> ("isDetachedFromSolidMesh");

  // We use the color map scheme because we can mark a face to be rupture ready from a partition where the face is a ghost.

  // process nodes on the interior
  {
    ModifiedObjectLists modifiedObjects;
    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
    {
      if( layersFromDomainBoundary[a]>1 &&
          (isSeparable[a])&&
          isNodeGhost[a]<0 &&
          nodeManager.m_toElementsRelation[a].size()>1 &&
          isDetachedFromSolidMesh[a] == 0 )
      {
        int splitMode = CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac);
        if (splitMode == 3)
        {
          ProcessKink(a, nodeManager, faceManager, modifiedObjects) ;
        }
        else if ( (splitMode == 0 || splitMode == 1) && nodeManager.m_toElementsRelation[a].size() >= 2)
        {
          ProcessTip(a, nodeManager, faceManager, splitMode, modifiedObjects) > 0;
        }
      }
    }
  }

  for( int color=0 ; color<partition.NumColor() ; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( partition.Color() == color )
    {

      // process "near-boundary" nodes
      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
      {
        if( layersFromDomainBoundary[a]<=1 &&
            (isSeparable[a]) &&
            isNodeGhost[a]<0 &&
            nodeManager.m_toElementsRelation[a].size()>1 &&
            isDetachedFromSolidMesh[a] == 0 )

        {
          int splitMode = CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac);
          if (splitMode == 3)
          {
            ProcessKink(a, nodeManager, faceManager, modifiedObjects) ;
          }
          else if ( (splitMode == 0 || splitMode == 1) && nodeManager.m_toElementsRelation[a].size() >= 2)
          {
            ProcessTip(a, nodeManager, faceManager, splitMode, modifiedObjects) > 0;
          }
        }
      }

    }
    partition.ModifyGhostsAndNeighborLists( modifiedObjects );
  } // color loop

}

//bool Fractunator2D::CheckAndSplitNode( const localIndex nodeID,
//                          NodeManagerT& nodeManager,
//                          EdgeManagerT& edgeManager,
//                          FaceManagerT& faceManager,
//                          ExternalFaceManagerT& externalFaceManager,
//                          ElementManagerT& elementManager,
//                          ModifiedObjectLists& modifiedObjects,
//                          const bool prefrac)
//{
//
//  // 0  We first check whether  node is a singular point (which should split regardless fracture criterion) or not.  If it is, CheckSplitability will returen 2 and we split it.
//  // 1. If a node is not a singular point, check whether it is a tip or not.  If it is a tip, calculate the stress intensity factors.
//  // 2. If it is a tip and SIF passes threshold, we first heal all faces connected to this node. We we pick the face with the greatest stress and mark its rupture state.
//  //      Even though UpdateRuptureStates already marked the rupture state, we don't use it.  We only use UpdateRuptureStates to calculate stressonface.  We will reevaluate the rupture state.
//  // 3. We then process this node and it should spit through the face that we have just marked.  If not, something is wrong.
//
//
//  int splitMode;
//  if (m_failCriterion==1)
//    {
//    splitMode = CheckNodeSplitability( nodeID, nodeManager, faceManager, edgeManager, prefrac);
//    }
//  bool isSplit = false;
//
//  if (splitMode == 2 || m_failCriterion == 0)
//  {
//    // In the former, this node is geometrically singular because of the splitting of a neighbor node.  We have to split it.  No need to check fracturing criterion.
//
//    isSplit= EvaluateAndSplitNode( nodeID, nodeManager, edgeManager, faceManager, externalFaceManager, elementManager, modifiedObjects ) ;
//  }
//  else if (splitMode == 3)
//  {
//    if (ProcessKink(nodeID, nodeManager, faceManager) > 0)
//    {
//      isSplit = EvaluateAndSplitNode( nodeID, nodeManager, edgeManager, faceManager, externalFaceManager, elementManager, modifiedObjects ) ;
//    }
//  }
//  else if ( (splitMode == 0 || splitMode == 1) && nodeManager.m_toElementsRelation[nodeID].size() >= 2)
//  {
//    if (ProcessTip(nodeID, nodeManager, faceManager, splitMode) > 0)
//    {
//      isSplit = EvaluateAndSplitNode( nodeID, nodeManager, edgeManager, faceManager, externalFaceManager, elementManager, modifiedObjects ) ;
//    }
//  }
//
//  return(isSplit);
//
//
//}




bool Fractunator2D::EvaluateAndSplitNode( const localIndex nodeID,
                                 NodeManagerT& nodeManager,
                                 EdgeManagerT& edgeManager,
                                 FaceManagerT& faceManager,
                                 ExternalFaceManagerT& externalFaceManager,
                                 ElementManagerT& elementManager,
                                 ModifiedObjectLists& modifiedObjects)
{

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );


  iArray1d& isDetachedFromSolidMesh = nodeManager.GetFieldData<int>("isDetachedFromSolidMesh");
  rArray1d& kinkAngle = nodeManager.GetFieldData<realT>("kinkAngle");
  rArray1d& SIF_I = nodeManager.GetFieldData<realT>("SIF_I");
  rArray1d& SIF_II = nodeManager.GetFieldData<realT>("SIF_II");

  lSet nodesForRenewKinkAngle;

  bool isSplit = false;

  if (isDetachedFromSolidMesh[nodeID] == 0)
  {

    iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");

    iArray1d* flowEdgeType = edgeManager.GetFieldDataPointer<int>("flowEdgeType");
    iArray1d* flowFaceType = faceManager.GetFieldDataPointer<int>("flowFaceType");
    rArray1d& stressNOnFace = faceManager.GetFieldData<realT>("stressNOnFace");
    Array1dT<R1Tensor>& stressTOnFace = faceManager.GetFieldData<R1Tensor>("stressTOnFace");

    Array1dT<lSet>& edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");


    //  nodeManager.AddKeylessDataField<int>( "nRupturedFace",true, true );
    //  nodeManager.AddKeylessDataField<int>( "nOpenFace",true, true );

    localIndex nRupturedFaces = 0, nExternalFaces = 0;
    lArray1d facesInvolved;

    for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
        iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
    {
      if (ruptureState[*iface] == 1)
      {
        nRupturedFaces++;
        facesInvolved.push_back(*iface);
      }
    }



    //We don't want to combine these loops because we want ruptured faces to be ahead of open faces in the list.
    for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
        iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
    {
      if (faceManager.m_isExternal[*iface] == 1) //(faceManager.m_toElementsRelation[*iface].size() == 1)
      {
        nExternalFaces++;
        facesInvolved.push_back(*iface);
      }
    }

    if ((nRupturedFaces >= 1 and nRupturedFaces + nExternalFaces >=2) or nExternalFaces >=4)
    {

      // Now we know that this node is going to be split.
      // To be safe (and lazy), we assume all the elements connected to this node are to be modified, and all faces, edges, and nodes belonging to these elements are to be modified.
      for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
          k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
      {
        localIndex iEle = k->second;
        ElementRegionT& elementRegion = *(k->first);

        modifiedObjects.modifiedElements[elementRegion.m_regionName].insert(iEle);

        for (localIndex j = 0; j < elementRegion.m_toNodesRelation.Dimension(1); ++j)
        {
          modifiedObjects.modifiedNodes.insert(elementRegion.m_toNodesRelation[iEle][j]);
          modifiedObjects.modifiedEdges.insert(elementRegion.m_toNodesRelation[iEle][j] );
        }

        for (localIndex j = 0; j < elementRegion.m_toFacesRelation.Dimension(1); ++j)
        {
          modifiedObjects.modifiedFaces.insert(elementRegion.m_toFacesRelation[iEle][j]);
        }

      }

      // We have to do this because in the case of hanging node (connected to 4 external faces), some faces are not connected to elements.
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
      {
        modifiedObjects.modifiedFaces.insert(*iface);
      }


      // To figure out along which two faces the node should split.
      if (nRupturedFaces == 2 && stressNOnFace[facesInvolved[0]] < stressNOnFace[facesInvolved[1]])
      {
        // We need to put the two faces with greatest tensile stresses at the beginning of the list.
        localIndex temp = facesInvolved[1];
        facesInvolved[1] = facesInvolved[0];
        facesInvolved[0] = temp;
      }
      else if ( nRupturedFaces > 2)
      {
        localIndex temp;
        for (localIndex i = 0; i < 2; ++i)
        {
          for (localIndex j = nRupturedFaces - 1; j > i ; --j)
          {
            if ( stressNOnFace[facesInvolved[j]] > stressNOnFace[facesInvolved[j-1]])
            {
              temp = facesInvolved[j-1];
              facesInvolved[j-1] = facesInvolved[j];
              facesInvolved[j] = temp;
            }
          }
        }
      }


      int nFace2Split = 0;
      localIndex dividingFaces[2];

      isSplit = true;
      localIndex daughterEdges[2];
      localIndex daughterNodes[2];

      SIF_I[nodeID] = 0.0;
      SIF_II[nodeID] = 0.0;
      nodeManager.SplitObject( nodeID, rank, daughterNodes);
      edgeManager.SplitObject( nodeID, rank, daughterEdges);


      if (m_verbose)
      {
        std::cout << "node " << nodeID << " split into " << daughterNodes[0]<< " and " << daughterNodes[1] << std::endl;
      }

      modifiedObjects.modifiedNodes.insert(nodeID);
      modifiedObjects.modifiedEdges.insert(nodeID);
      modifiedObjects.newNodes.insert(daughterNodes[0]);
      modifiedObjects.newNodes.insert(daughterNodes[1]);
      modifiedObjects.newEdges.insert(daughterEdges[0]);
      modifiedObjects.newEdges.insert(daughterEdges[1]);

      nodesForRenewKinkAngle.insert(daughterNodes[0]);
      nodesForRenewKinkAngle.insert(daughterNodes[1]);
      isDetachedFromSolidMesh[nodeID] = 1;
      kinkAngle[nodeID] = -1;

      if ( flowEdgeType != NULL) (*flowEdgeType)[nodeID] = 0;

      edgeManager.m_toNodesRelation(daughterNodes[0],0) = daughterNodes[0];
      edgeManager.m_toNodesRelation(daughterNodes[1],0) = daughterNodes[1];

      nodeManager.m_nodeToEdgeMap[daughterNodes[0]].insert(daughterNodes[0]);
      nodeManager.m_nodeToEdgeMap[daughterNodes[1]].insert(daughterNodes[1]);

      //      (*nodeManager.m_refposition)[daughterNodes[0]][2] += ((realT(rand()) / RAND_MAX) ) * 1.0;
      //      (*nodeManager.m_refposition)[daughterNodes[1]][2] -= ((realT(rand()) / RAND_MAX) ) * 1.0;

      if (nodeManager.m_mass != NULL)
      {
        (*nodeManager.m_mass)[daughterNodes[0]] = 0.5 * (*nodeManager.m_mass)[nodeID] ;
        (*nodeManager.m_mass)[daughterNodes[1]] = 0.5 * (*nodeManager.m_mass)[nodeID] ;
        (*nodeManager.m_mass)[nodeID] = 1e-99;
      }


      if (nExternalFaces >=4)
        // In this case we just have four external faces connected to this node.
        // The two faces that we are going to pick cannot be sisters and there must be at least one element between them.
      {
        nFace2Split = 0;
        //lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;

        dividingFaces [0] = facesInvolved[nRupturedFaces];  // This should be the first external face in the list

        localIndex i = nRupturedFaces;

        unsigned int numEleOnOneSide = 0;

        while (numEleOnOneSide == 0)
        {
          i++;

          if (i >= nRupturedFaces + nExternalFaces)
            throw GPException("The node has four external faces but we cannot find a viable path to split it.");
          dividingFaces [1] = facesInvolved[i];

          R1Tensor vecWings[2], vecWingNorm;

          int iDivType = MakeDividingVectors ( nodeID, vecWingNorm, vecWings, dividingFaces, nodeManager, faceManager);

          for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
              k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
          {
            int mySide = MySideAlongFracture(nodeID, *k, iDivType, vecWingNorm, vecWings,
                                             nodeManager, elementManager);

            numEleOnOneSide += mySide;
          }

          if (numEleOnOneSide == nodeManager.m_toElementsRelation[nodeID].size()) numEleOnOneSide = 0;
        }
      }
      else if ( nRupturedFaces >= 1 && nExternalFaces >= 1)
        //In this case we split along one of the external faces and the first ruptured face (after sorting, this face should have the greatest normal stress) .
      {
        nFace2Split = 1;
        dividingFaces [0] = facesInvolved [0];
        dividingFaces [1] = facesInvolved [nRupturedFaces + nExternalFaces -1 ];

      }
      else //if ( nRupturedFaces >= 1)
        // In this case we just pick any two ruptured faces to split.  The other ruptured faces will be split in a later iteration.
      {
        nFace2Split = 2;
        dividingFaces [0] = facesInvolved [0];
        dividingFaces [1] = facesInvolved [1];

      }


      int iDivType;
      // iDivType = 1 if the two wings of the dividing faces from a straight line; =2 if they form an angle.
      R1Tensor vecWings[2], vecWingNorm;

      iDivType = MakeDividingVectors ( nodeID, vecWingNorm, vecWings, dividingFaces, nodeManager, faceManager);


      // Now we split the ruptured faces.
      // Here we only handles the connectivity between faces and elements.  We will deal with nodes later.

      rArray1d* delta0N = faceManager.GetFieldDataPointer<realT>("delta0N");
      Array1dT<R1Tensor>* stressShear0 = faceManager.GetFieldDataPointer<R1Tensor>("stressShear0");

      for (int i = 0; i < nFace2Split; ++i)
      {
        localIndex iFace = facesInvolved[i];
        localIndex daughterFaces[2];
        faceManager.SplitObject( iFace, rank, daughterFaces);

        modifiedObjects.modifiedFaces.insert(iFace);
        modifiedObjects.newFaces.insert(daughterFaces[0]);
        modifiedObjects.newFaces.insert(daughterFaces[1]);

        //TODO: add logic for updating contact manager, external face manager, neighborlist, etc.

        ruptureState[iFace] = 2;
        ruptureState[daughterFaces[0]] = 3;
        ruptureState[daughterFaces[1]] = 3;

        if (delta0N != NULL)
        {
          //realT aFace = faceManager.SurfaceArea(nodeManager, iFace, true);
          if (stressNOnFace[iFace] < 0)
          {
            // We pass the contact stiffness in from the solver via field dleta0 itself.
            (*delta0N)[iFace] = -stressNOnFace[iFace] /(*delta0N)[iFace];
            (*stressShear0)[iFace] = stressTOnFace[iFace];
          }
          else
          {
            (*delta0N)[iFace] = 0.0;
            (*stressShear0)[iFace] = 0.0;
          }
        }

        localIndex theOtherNode;
        theOtherNode = faceManager.m_toNodesRelation[iFace][0] + faceManager.m_toNodesRelation[iFace][1] - nodeID;
        nodesForRenewKinkAngle.insert(theOtherNode);


        if ( flowFaceType != NULL ) (*flowFaceType)[iFace] = 0;

        for( lArray1d::iterator ie = faceManager.m_toEdgesRelation[iFace].begin() ;
            ie!=faceManager.m_toEdgesRelation[iFace].end() ; ++ie )
        {
          localIndex edgeID = *ie;
          while (edgeManager.m_parentIndex[edgeID] < edgeManager.DataLengths())
          {
            edgeID = edgeManager.m_parentIndex[edgeID];
            modifiedObjects.modifiedEdges.insert(edgeID);
          }

          edgeToFlowFaces[edgeID].insert(iFace);
        }



        faceManager.m_isExternal[daughterFaces[0]] = 1;
        faceManager.m_isExternal[daughterFaces[1]] = 1;

        faceManager.m_toNodesRelation[daughterFaces[0]] = faceManager.m_toNodesRelation[iFace];
        faceManager.m_toNodesRelation[daughterFaces[1]] = faceManager.m_toNodesRelation[iFace];

        faceManager.m_toEdgesRelation[daughterFaces[0]] = faceManager.m_toEdgesRelation[iFace];
        faceManager.m_toEdgesRelation[daughterFaces[1]] = faceManager.m_toEdgesRelation[iFace];

        externalFaceManager.SplitFace(daughterFaces[0], daughterFaces[1], nodeManager);

        for (localIndex a = 0; a < 2; ++a)
        {
          localIndex jFace = daughterFaces[a];
          for( lArray1d::iterator b = faceManager.m_toNodesRelation[jFace].begin() ;
              b!=faceManager.m_toNodesRelation[jFace].end() ; ++b )
          {
            if (*b != nodeID)
            {
              nodeManager.m_nodeToFaceMap[*b].insert(jFace);
              edgeManager.m_toFacesRelation[*b].insert(jFace);
            }
          }
        }

        // Connect new faces to the correct new nodes and new edges.
        for ( localIndex a = 0; a < 2; ++a)
        {
          for( lArray1d::iterator b = faceManager.m_toNodesRelation[daughterFaces[a]].begin() ;
              b!=faceManager.m_toNodesRelation[daughterFaces[a]].end() ; ++b )
          {
            if (*b == nodeID) *b = daughterNodes[a];
          }

          for( lArray1d::iterator b = faceManager.m_toEdgesRelation[daughterFaces[a]].begin() ;
              b!=faceManager.m_toEdgesRelation[daughterFaces[a]].end() ; ++b )
          {
            if (*b == nodeID) *b = daughterEdges[a];
          }

        }

        // Connect the dead face to the correct nodes.  If the node to split is connected to a dead face, then the dead face should attache to the node's up-most ancestor.
        for( lArray1d::iterator b = faceManager.m_toNodesRelation[iFace].begin() ;
            b!=faceManager.m_toNodesRelation[iFace].end() ; ++b )
        {
          while (nodeManager.m_parentIndex[*b] < nodeManager.DataLengths())
          {
            nodeManager.m_nodeToFaceMap[*b].erase(iFace);
            *b = nodeManager.m_parentIndex[*b];
            nodeManager.m_nodeToFaceMap[*b].insert(iFace);
          }
        }

        for( lArray1d::iterator b = faceManager.m_toEdgesRelation[iFace].begin() ;
            b!=faceManager.m_toEdgesRelation[iFace].end() ; ++b )
        {
          if (*b==nodeID)
          {
            while (edgeManager.m_parentIndex[*b] < edgeManager.DataLengths())
            {
              edgeManager.m_toFacesRelation[*b].erase(iFace);
              *b = edgeManager.m_parentIndex[*b];
              edgeManager.m_toFacesRelation[*b].insert(iFace);
            }

          }
        }



        localIndex iSide;
        for (Array1dT< std::pair< ElementRegionT*, localIndex > >::iterator iter = faceManager.m_toElementsRelation[iFace].begin();
            iter != faceManager.m_toElementsRelation[iFace].end(); ++iter )
        {
          ElementRegionT* elementRegion = iter->first;
          localIndex iEle= iter->second;

          iSide = MySideAlongFracture (nodeID, *iter , iDivType,
                                       vecWingNorm, vecWings,
                                       nodeManager, elementManager );

          localIndex* facelist = elementRegion->m_toFacesRelation[iEle];
          for (localIndex a=0 ; a<elementRegion->m_toFacesRelation.Dimension(1) ; ++a )
          {
            localIndex jFace = facelist[a];
            if (jFace == iFace )
            {
              facelist[a] = faceManager.m_childIndices[iFace][iSide];
              faceManager.m_toElementsRelation[facelist[a]].resize(1);
              faceManager.m_toElementsRelation[facelist[a]][0].first = elementRegion;
              faceManager.m_toElementsRelation[facelist[a]][0].second = iEle;
            }

          }

        }

        for (Array1dT< std::pair< ElementRegionT*, localIndex > >::iterator iter = faceManager.m_toElementsRelation[iFace].end();
            iter != faceManager.m_toElementsRelation[iFace].begin(); --iter)
        {
          faceManager.m_toElementsRelation[iFace].erase(iter);
        }

        //faceManager.m_toElementsRelation[iFace].resize(0);

        faceManager.SortFaceNodes( nodeManager, daughterFaces[0] );
        faceManager.SortFaceNodes( nodeManager, daughterFaces[1] );

        // We need to figure out the sign of initial stress based on which child face is along the original positive normal direction
        {
          R1Tensor n0, n1;
          n0 = faceManager.FaceNormal(nodeManager, iFace, true);
          n1 = faceManager.FaceNormal(nodeManager, daughterFaces[0], true);
          if (Dot(n0, n1) < 0.0 && stressShear0 != NULL)
          {
            (*stressShear0)[iFace] *= -1;
          }

        }
      }


      // Loop through all elements connected to the node to be split.
      for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
          k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
      {
        localIndex iEle = k->second;
        ElementRegionT* elementRegion = k->first;

        int mySide = MySideAlongFracture(nodeID, *k, iDivType, vecWingNorm, vecWings,
                                         nodeManager, elementManager);

        nodeManager.m_toElementsRelation[daughterNodes[mySide]].insert(*k);

        // Correct the element to node map.
        for (localIndex iN = 0; iN < elementRegion->m_toNodesRelation.Dimension(1); ++iN)
        {
          if ( elementRegion->m_toNodesRelation(iEle,iN) == nodeID)
          {
            elementRegion->m_toNodesRelation(iEle,iN) = daughterNodes[mySide];
          }
        }

        // Correct the face to node map.
        for (localIndex iF = 0; iF < elementRegion->m_toFacesRelation.Dimension(1); ++iF)
        {
          localIndex faceID = elementRegion->m_toFacesRelation(iEle,iF);
          for( lArray1d::iterator b = faceManager.m_toNodesRelation[faceID].begin() ;
              b!=faceManager.m_toNodesRelation[faceID].end() ; ++b )
          {
            if (*b == nodeID)
            {
              *b = daughterNodes[mySide];
              nodeManager.m_nodeToFaceMap[nodeID].erase(faceID);
            }
          }
          // Correct the face to edge map.
          for( lArray1d::iterator b = faceManager.m_toEdgesRelation[faceID].begin() ;
              b!=faceManager.m_toEdgesRelation[faceID].end() ; ++b )
          {
            if (*b == nodeID)
            {
              *b = daughterEdges[mySide];
              edgeManager.m_toFacesRelation[nodeID].erase(faceID);
            }
          }
        }

      }
      nodeManager.m_toElementsRelation[nodeID].clear();


      // Build node to face map for new nodes

      for (localIndex iD = 0; iD < 2; ++iD)
      {
        localIndex iNd = daughterNodes[iD];
        localIndex iEdge = daughterEdges[iD];

        nodeManager.m_nodeToFaceMap[iNd].clear();
        edgeManager.m_toFacesRelation[iEdge].clear();


        // Loop through all elements connected to this new node.
        for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[iNd].begin() ;
            k!=nodeManager.m_toElementsRelation[iNd].end() ; ++k )
        {
          localIndex iEle = k->second;
          ElementRegionT* elementRegion = k->first;

          localIndex* facelist = elementRegion->m_toFacesRelation[iEle];
          for (localIndex a=0 ; a<elementRegion->m_toFacesRelation.Dimension(1) ; ++a )
          {
            localIndex faceID = facelist[a];

            bool faceRelatedToThisNode = false;
            for( lArray1d::const_iterator b=faceManager.m_toNodesRelation[faceID].begin() ;
                b!=faceManager.m_toNodesRelation[faceID].end() ; ++b )
            {
              if (*b == iNd) faceRelatedToThisNode = true;
            }

            if (faceRelatedToThisNode && nodeManager.m_nodeToFaceMap[iNd].count(faceID) == 0)
            {
              nodeManager.m_nodeToFaceMap[iNd].insert(faceID);
            }
            if (faceRelatedToThisNode && edgeManager.m_toFacesRelation[iEdge].count(faceID) == 0)
            {
              edgeManager.m_toFacesRelation[iEdge].insert(faceID);
              //              if ( flowFaceType != NULL )
              //              {
              //                if ( (*flowFaceType)[faceID] > -1)
              //                {
              //                  edgeToFlowFaces[iEdge].insert(faceID);
              //                }
              //              }

            }
          }
        }
      }

      // This is causing problems in parallel
      //      if (!prefrac) // This is to make sure one node is only split once in one time step.
      //      {
      //        for (localIndex i = 0; i < 2; ++i)
      //        {
      //          for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[daughterNodes[i]].begin() ;
      //              iface!=nodeManager.m_nodeToFaceMap[daughterNodes[i]].end() ; ++iface )
      //          {
      //            if (ruptureState[*iface] == 1)
      //
      //            {
      //              std::cout << "resetting rupture state of face " << *iface << "at node " << daughterNodes[i] << std::endl;
      //              ruptureState[*iface] = 0;
      //            }
      //          }
      //
      //        }
      //
      //      }

    }


  } // if this is not a dead node.



  if (isSplit)
  {
    for( lSet::const_iterator i=nodesForRenewKinkAngle.begin() ; i!=nodesForRenewKinkAngle.end() ; ++i )
    {
      kinkAngle[*i] = CalculateKinkAngle(*i, nodeManager, faceManager);
    }
  }

  if (isSplit && false)
  {
    std::ofstream outfile_;
    string fileName;
    fileName = "T";
    fileName.append(static_cast<ostringstream*>( &(ostringstream() << nodeID) )->str());
    fileName.append(".txt");
    outfile_.open(fileName.c_str());
    WriteTopology (outfile_, nodeManager, elementManager, faceManager, edgeManager);
  }

  if (isSplit)
  {
    CorrectSplitNodalMass(nodeManager,
                          nodeManager.m_childIndices[nodeID][0],
                          nodeManager.m_childIndices[nodeID][1]);
  }
  return isSplit;

}





int Fractunator2D::MakeDividingVectors (localIndex nodeID,
                                        R1Tensor& vecWingNorm,
                                        R1Tensor vecWings[],
                                        localIndex dividingFaces[],
                                        NodeManagerT& nodeManager,
                                        FaceManagerT& faceManager)
{
  int iDivType;
  faceManager.FaceVector2D( nodeManager, dividingFaces[0], nodeID, vecWings[0]);
  faceManager.FaceVector2D( nodeManager, dividingFaces[1], nodeID, vecWings[1]);

  realT a1 = vecWings[0].L2_Norm();
  realT a2 = vecWings[1].L2_Norm();
  realT a0 = Dot(vecWings[0], vecWings[1]);
  if (a0 < -0.99999 * a1 * a2)
  {
    iDivType = 1;
  }
  else
  {
    iDivType = 2;
  }

  if (iDivType == 1)
  {
    // In this scenario, we only need the norm vector of the line constituted by the two wings.
    // We don't care about the length of this vector.  It can point to either side of the line.
    vecWingNorm[0] = -vecWings[0][1];
    vecWingNorm[1] = vecWings[0][0];
    vecWingNorm[2] = 0.0;
  }
  else
  {
    //We flip the two wing vectors if the angle from the first one to the second one is greater than 180 degrees.
    if ( vecWings[0][0] * vecWings[1][1] - vecWings[0][1] * vecWings[1][0] < 0.0 )
    {
      R1Tensor vec = vecWings[0];
      vecWings[0] = vecWings[1];
      vecWings[1] = vec;
    }

    vecWings[0].Normalize();
    vecWings[1].Normalize();

  }

  return iDivType;

}

int Fractunator2D::MySideAlongFracture (localIndex nodeID,
                                        const std::pair< ElementRegionT*, localIndex >& elem,
                                        int iDivType,
                                        R1Tensor vecWingNorm,
                                        R1Tensor vecWings[],
                                        NodeManagerT& nodeManager,
                                        ElementManagerT& elementManager )
{
  int mySide;
  R1Tensor eleCenter(0.0, 0.0, 0.0);
  const ElementRegionT& elemRegion = *(elem.first);
  localIndex eleID = elem.second;

  for (localIndex i = 0; i != elemRegion.m_toNodesRelation.Dimension(1); ++i)
  {
    localIndex iNd = elemRegion.m_toNodesRelation[eleID][i];
    eleCenter = eleCenter + (*nodeManager.m_refposition)[iNd];

  }

  eleCenter = eleCenter / elemRegion.m_toNodesRelation.Dimension(1);
  eleCenter = eleCenter -  (*nodeManager.m_refposition)[nodeID];
  eleCenter.Normalize();

  if (iDivType == 1)
  {
    if (Dot(eleCenter, vecWingNorm) > 0)
    {
      mySide = 0;
    }
    else
    {
      mySide = 1;
    }
  }
  else
  {
    if ( vecWings[0][0] * eleCenter [1] - vecWings[0][1] * eleCenter [0] > 0.0 &&
        Dot(vecWings[0], eleCenter) > Dot(vecWings[0], vecWings[1]))
    {
      mySide = 0;
    }
    else
    {
      mySide = 1;
    }
  }
  return mySide;
}


void Fractunator2D::FindClosestNode (NodeManagerT& nodeManager,
                                     R1Tensor& xPt,
                                     localIndex& iNdMinDist )
{
  realT minDist, dist;
  R1Tensor v;
  v = xPt;
  v -= (*nodeManager.m_refposition)[0];
  minDist = v.L2_Norm();
  iNdMinDist = 0;
  for (localIndex iNd = 1; iNd < (nodeManager.DataLengths()) ; ++iNd)
  {
    v = xPt;
    v -= (*nodeManager.m_refposition)[iNd];
    dist = v.L2_Norm();
    if ( dist < minDist )
    {
      minDist = dist;
      iNdMinDist = iNd;
    }
  }
}

void Fractunator2D::MarkPreexistingFractures (FaceManagerT& faceManager,
                                              NodeManagerT& nodeManager,
                                              SpatialPartition& partition,
                                              ModifiedObjectLists& modifiedObjects)
{

  partition.UpdatePartitionBoundingBox(nodeManager);
  for (localIndex iPF = 0; iPF < m_x1_PreFrac.size(); ++iPF)
  {
    bool lineInPartition = false;
    int nPtInPartition, nIntersection;
    R1Tensor xEnd[2], xInterEdge[8];

    nPtInPartition = 0;
    if ( partition.IsCoordInPartitionClosed(m_x0_PreFrac[iPF]) )
    {
      xInterEdge[nPtInPartition] = m_x0_PreFrac[iPF];
      nPtInPartition = 1;
      lineInPartition = true;
    }
    if ( partition.IsCoordInPartitionClosed(m_x1_PreFrac[iPF]) )
    {
      xInterEdge[nPtInPartition] = m_x1_PreFrac[iPF];
      nPtInPartition += 1;
      lineInPartition = true;
    }


    R1TensorT<2> x0PartitionEdge[4], x1PartitionEdge[4], x0Line, x1Line;
    R1Tensor xMin, xMax;

    partition.getPartitionGeometricalBoundary(xMin, xMax);


    x0PartitionEdge[0][0] = xMin[0];
    x0PartitionEdge[0][1] = xMin[1];
    x1PartitionEdge[0][0] = xMax[0];
    x1PartitionEdge[0][1] = xMin[1];

    x0PartitionEdge[1][0] = xMax[0];
    x0PartitionEdge[1][1] = xMin[1];
    x1PartitionEdge[1][0] = xMax[0];
    x1PartitionEdge[1][1] = xMax[1];

    x0PartitionEdge[2][0] = xMax[0];
    x0PartitionEdge[2][1] = xMax[1];
    x1PartitionEdge[2][0] = xMin[0];
    x1PartitionEdge[2][1] = xMax[1];

    x0PartitionEdge[3][0] = xMin[0];
    x0PartitionEdge[3][1] = xMax[1];
    x1PartitionEdge[3][0] = xMin[0];
    x1PartitionEdge[3][1] = xMin[1];

    x0Line[0] = m_x0_PreFrac[iPF][0];
    x0Line[1] = m_x0_PreFrac[iPF][1];
    x1Line[0] = m_x1_PreFrac[iPF][0];
    x1Line[1] = m_x1_PreFrac[iPF][1];

    R1TensorT<2> xInter[2];
    nIntersection = 0;

    for (localIndex i = 0; i < 4; ++i)
    {
      localIndex interCase = GeometryUtilities::LineIntersection(x0PartitionEdge[i],
                                                                 x1PartitionEdge[i],
                                                                 x0Line,
                                                                 x1Line,
                                                                 xInter[0],
                                                                 xInter[1],
                                                                 -1.0e-6);
      if (interCase ==0 || interCase == 1)
      {
        lineInPartition = true;
        xInterEdge[nIntersection + nPtInPartition][0] = xInter[0][0];
        xInterEdge[nIntersection + nPtInPartition][1] = xInter[0][1];
        nIntersection++;
      }
      else if (interCase == 2)
      {
        lineInPartition = true;
        xInterEdge[nIntersection + nPtInPartition][0] = xInter[0][0];
        xInterEdge[nIntersection + nPtInPartition][1] = xInter[0][1];
        nIntersection++;
        xInterEdge[nIntersection + nPtInPartition][0] = xInter[1][0];
        xInterEdge[nIntersection + nPtInPartition][1] = xInter[1][1];
        nIntersection++;
      }
    }

    if (nPtInPartition == 2)
    {
      xEnd[0] = m_x0_PreFrac[iPF];
      xEnd[1] = m_x1_PreFrac[iPF];
    }
    else if (nPtInPartition == 1 && nIntersection ==1)
    {
      xEnd[0] = xInterEdge[0];
      xEnd[1] = xInterEdge[1];
    }
    else if (nPtInPartition + nIntersection >= 2)
    {
      realT maxL = 0.0;
      localIndex ii,jj;
      ii = 0;
      jj = 1;
      for (int i = 0; i < nPtInPartition + nIntersection -1; ++i)
      {
        for (int j=i+1; j < nPtInPartition + nIntersection; ++j)
        {
          if ( (xInterEdge[i]-xInterEdge[j]).Normalize() > maxL)
          {
            maxL = (xInterEdge[i]-xInterEdge[j]).Normalize();
            ii = i;
            jj = j;
          }
        }
      }
      xEnd[0] = xInterEdge[ii];
      xEnd[1] = xInterEdge[jj];
    }
    else if ( nPtInPartition + nIntersection <= 1 )
    {
      //throw GPException("Error in finding intersections between preexisting fractures and partitions");
      // Due to the error in checking line intersection, nPtInPartion + nIntersection could be one.
      // We just ignore this case.
      lineInPartition = false;
    }

    if (lineInPartition)
    {
      //      xEnd[0] += xEnd[1];
      //      xEnd[0] /= 2.0;

      //      xEnd[1] = m_x0_PreFrac[iPF];
      //      TrackAlongAFractures(faceManager,
      //                           nodeManager,
      //                           xEnd[0],
      //                           xEnd[1],
      //                           modifiedObjects);

      //      TrackAlongAFractures(faceManager,
      //                           nodeManager,
      //                           xEnd[1],
      //                           xEnd[0],
      //                           modifiedObjects);

      TrackAlongShortestPath(faceManager,
                             nodeManager,
                             xEnd[0],
                             xEnd[1],
                             modifiedObjects);

    }

  }

}

void Fractunator2D::TrackAlongAFractures (FaceManagerT& faceManager,
                                          NodeManagerT& nodeManager,
                                          R1Tensor x0,
                                          R1Tensor x1,
                                          ModifiedObjectLists& modifiedObjects)
{
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  localIndex iClosestNd;
  R1Tensor xFracFront, R1;
  realT r1, r2, r3, nDist;
  FindClosestNode( nodeManager, x0, iClosestNd);


  GeometryUtilities::ProjectPointToLineSegment( x0,
                                                x1,
                                                (*nodeManager.m_refposition)[iClosestNd],
                                                r1, r2, r3,
                                                xFracFront);
  localIndex iNdPFFront = iClosestNd;
  int flagDone = 0;
  localIndex iFacePicked, iNdPicked;
  realT lPrjFaceOnPF;
  lSet iNdPickedEver;


  while (flagDone == 0)
  {
    //A vector pointing from the current fracture front to the expected end of the preexisting fracture.
    R1Tensor xPFVec = x1;
    R1Tensor xFaceVec;

    xPFVec -= (*nodeManager.m_refposition)[iNdPFFront];
    // The do while loop continues only if a new eligible node can be found along the PF.
    flagDone=1;
    realT lPF = xPFVec.L2_Norm();
    realT minDist = std::numeric_limits<double>::max();
    iFacePicked = 0;

    iNdPicked = 0;

    if (lPF > 0 )
    {
      for( lSet::const_iterator iFace=nodeManager.m_nodeToFaceMap[iNdPFFront].begin() ;
          iFace!=nodeManager.m_nodeToFaceMap[iNdPFFront].end() ; ++iFace )
      {
        faceManager.FaceVector2D( nodeManager, *iFace, iNdPFFront, xFaceVec);
        // Length of the projection of this face (edge) on the preexisting fracture;
        lPrjFaceOnPF = xFaceVec * xPFVec / lPF;

        if ( lPrjFaceOnPF > 0.0 && faceManager.m_isExternal[*iFace] == 0)
        {
          flagDone = 0; // Now we know there is at least one qualified edge.  We need to choose one.

          localIndex iTheOtherNd = faceManager.m_toNodesRelation[*iFace][0] +
              faceManager.m_toNodesRelation[*iFace][1] - iNdPFFront;

          GeometryUtilities::ProjectPointToLineSegment( (*nodeManager.m_refposition)[iNdPFFront],
                                                        x1,
                                                        (*nodeManager.m_refposition)[iTheOtherNd],
                                                        nDist, r2, r3,
                                                        R1);
          if ( fabs(nDist) / lPrjFaceOnPF < minDist)
          {
            minDist = fabs(nDist) / lPrjFaceOnPF;
            iFacePicked = *iFace;
            iNdPicked = iTheOtherNd;

          }
        }
      } //Now we have searched all the faces connected to iNdPFFront.

      if (iNdPickedEver.count(iNdPicked) > 0)  // Sometimes the track point will go back and forth along the same face.
      {
        flagDone = 1;
      }
      else
      {
        iNdPickedEver.insert(iNdPicked);
      }

      xFaceVec = (*nodeManager.m_refposition)[iNdPicked];
      xFaceVec -= (*nodeManager.m_refposition)[iNdPFFront];
      lPrjFaceOnPF = xFaceVec * xPFVec / lPF;

      if (lPrjFaceOnPF > 1.2 * lPF) flagDone =1;  //The remaining part of the PF is not worth extending;

      if (flagDone == 0)
      {
        GeometryUtilities::ProjectPointToLineSegment( xFracFront,
                                                      x1,
                                                      (*nodeManager.m_refposition)[iNdPicked],
                                                      nDist, r2, r3,
                                                      R1);
        //Now we move the fracture front to the new node identified and prepare for the next advancement.
        xFracFront = R1;
        iNdPFFront = iNdPicked;

        //We avoid the situation where a single element could be separated from all its neighbor elements.
        int isOrphan = CheckOrphanElement ( faceManager, iFacePicked);

        if ( ruptureState[iFacePicked] == 0 && isOrphan == 0 ) //&& isGhost[iFacePicked] < 0 )
        {
          ruptureState[iFacePicked] = 1;
          modifiedObjects.modifiedFaces.insert(iFacePicked);
        }

        xPFVec = x1;
        xPFVec -= xFracFront;
        lPF = xPFVec.L2_Norm();

        if (fabs(nDist) > 0.5 * lPF) flagDone=1;

      }
    }
  }
}
//Find and mark the shortest path from node0 to node1 with Dijkstra's algorithm

void Fractunator2D::TrackAlongShortestPath (FaceManagerT& faceManager,
                                            NodeManagerT& nodeManager,
                                            R1Tensor x0,
                                            R1Tensor x1,
                                            ModifiedObjectLists& modifiedObjects)
{
  const realT LARGEDISTANCE = 1.0e100;
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");

  rArray1d faceArea(faceManager.DataLengths());

  iArray1d nodeVisited(nodeManager.DataLengths());
  nodeVisited = 0;


  rArray1d distToNode0(nodeManager.DataLengths());
  distToNode0 = LARGEDISTANCE;

  for (localIndex i = 0; i != faceManager.DataLengths(); ++i)
  {
    R1Tensor faceCenter, prj;
    realT ndist1, ndist2, udist, seg;
    faceManager.FaceCenter(nodeManager, i, faceCenter);
    GeometryUtilities::ProjectPointToLineSegment(x0, x1, (*nodeManager.m_refposition)[faceManager.m_toNodesRelation[i][0]], ndist1, udist, seg, prj);
    GeometryUtilities::ProjectPointToLineSegment(x0, x1, (*nodeManager.m_refposition)[faceManager.m_toNodesRelation[i][1]], ndist2, udist, seg, prj);


    faceArea[i] = faceManager.SurfaceArea(nodeManager, i, true) + std::max(fabs(ndist1),fabs(ndist2));

  }

  localIndex node0, node1;
  FindClosestNode(nodeManager, x0, node0);
  FindClosestNode(nodeManager, x1, node1);

  localIndex currentNode = node0;
  distToNode0[node0] = 0.0;

  realT shortest = 0.0;
  while (nodeVisited[node1] != 1 && shortest != LARGEDISTANCE)
  {
    nodeVisited[currentNode] = 1;
    for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[currentNode].begin() ;
        iface!=nodeManager.m_nodeToFaceMap[currentNode].end() ; ++iface )
    {
      localIndex nodej = faceManager.m_toNodesRelation[*iface][0] + faceManager.m_toNodesRelation[*iface][1] - currentNode;
      distToNode0[nodej] = std::min(distToNode0[nodej], distToNode0[currentNode] + faceArea[*iface]);
    }

    localIndex nNdSearched = 0;

    for (localIndex iNd = 0; iNd < nodeManager.DataLengths(); ++iNd)
    {
      if (nodeVisited[iNd] == 0)
      {
        if (nNdSearched == 0)
        {
          currentNode = iNd;
          shortest = distToNode0[iNd];
          ++nNdSearched;
        }
        else if ( distToNode0[iNd] < shortest)
        {
          shortest = distToNode0[iNd];
          currentNode = iNd;
          ++nNdSearched;
        }
      }
    }
  }

  if (nodeVisited[node1] == 1) // We have find the shortest path from node0 to node1.  Now we track back and mark the faces.
  {
    currentNode = node1;
    while (currentNode != node0)
    {
      shortest = LARGEDISTANCE;
      int fracFace = -1;
      localIndex nextNode = currentNode;
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[currentNode].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[currentNode].end() ; ++iface )
      {
        localIndex nodej = faceManager.m_toNodesRelation[*iface][0] + faceManager.m_toNodesRelation[*iface][1] - currentNode;
        if (distToNode0[nodej] < shortest)
        {
          shortest = distToNode0[nodej];
          nextNode = nodej;
          fracFace = *iface;
        }
      }
      if (fracFace >= 0)
      {
        currentNode = nextNode;
        if (faceManager.m_toElementsRelation[fracFace].size() == 2)
        {
          ruptureState[fracFace] = 1;
          modifiedObjects.modifiedFaces.insert(fracFace);
        }
      }
    }
  }

}

void Fractunator2D::CreatePreexistingFractures (NodeManagerT& nodeManager,
                                                EdgeManagerT& edgeManager,
                                                FaceManagerT& faceManager,
                                                ExternalFaceManagerT& externalFaceManager,
                                                ElementManagerT& elementManager,
                                                SpatialPartition& partition,
                                                const bool prefrac )
{
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //  const iArray1d& layersFromDomainBoundary = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
  //  std::cout << "My color " << partition.Color() << " out of " << partition.NumColor() << std::endl;
  //  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");


  int nNdSplitThisIter = 100;
  //NodeManagerT* NodeManager = &nodeManager;
  //Array1dT<realT>& mass = NodeManager->GetFieldData<FieldInfo::mass> ();
  localIndex iIter=0;

  Array1dT<lSet> nodesToRupturedFaces;
  Array1dT<lSet> edgesToRupturedFaces;

  const iArray1d& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();

  if (rank==0) std::cout << "Preexisting fracture creation: iteration " << iIter << std::endl;

  for( int color=0 ; color<partition.NumColor() ; ++color )
  {

    ModifiedObjectLists modifiedObjects;

    if ( partition.Color() == color)
    {
      MarkPreexistingFractures(faceManager, nodeManager, partition, modifiedObjects); // We don't need to mark the fractures on subsequent iterations.  These iteratures are just to break nodes that are required to break by geometrical conditions.
    }
    partition.ModifyGhostsAndNeighborLists( modifiedObjects );

  }

  while (nNdSplitThisIter > 0)
  {
    nNdSplitThisIter = 0;

    //std::cout << rank<< "My color " << partition.Color() << std::endl;

    for( int color=0 ; color<partition.NumColor() ; ++color )
    {

      ModifiedObjectLists modifiedObjects;

      if ( partition.Color() == color)
      {


        for (localIndex iNd = 0; iNd < nodeManager.DataLengths(); ++iNd)
        {
          if (nodeManager.m_childIndices[iNd].size() == 0 &&   //A node can only be split once.
              isNodeGhost[iNd] < 0)

          {

            bool didSplit = EvaluateAndSplitNode(iNd,  nodeManager,
                                        edgeManager,
                                        faceManager,
                                        externalFaceManager,
                                        elementManager,
                                        modifiedObjects);


            if (didSplit)
            {
              nNdSplitThisIter +=1;
              std::cout << "Partition " << rank << " split node " << iNd << std::endl;
            }

          }
        }

      }

      partition.ModifyGhostsAndNeighborLists( modifiedObjects );

    }



    if ( nNdSplitThisIter > 0) std::cout << "Partition " << rank << " split " << nNdSplitThisIter << " nodes in this iteration." << std::endl;

    int myNSplit = nNdSplitThisIter;
    MPI_Allreduce(&myNSplit, &nNdSplitThisIter, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
    iIter++;
  }


}





bool Fractunator2D::FindFracturePlanes( const localIndex nodeID,
                                        const NodeManagerT& nodeManager,
                                        const EdgeManagerT& edgeManager,
                                        const FaceManagerT& faceManager,
                                        const Array1dT<lSet>& nodesToRupturedFaces,
                                        const Array1dT<lSet>& edgesToRupturedFaces,
                                        lSet& separationPathFaces,
                                        std::map<localIndex,int>& edgeLocations,
                                        std::map<localIndex,int>& faceLocations,
                                        std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{
  return (false);
}

void Fractunator2D::PerformFracture( const localIndex nodeID,
                                     NodeManagerT& nodeManager,
                                     EdgeManagerT& edgeManager,
                                     FaceManagerT& faceManager,
                                     ExternalFaceManagerT& externalFaceManager,
                                     ElementManagerT& elementManager,
                                     ModifiedObjectLists& modifiedObjects,
                                     Array1dT<lSet>& nodesToRupturedFaces ,
                                     Array1dT<lSet>& edgesToRupturedFaces ,
                                     const lSet& separationPathFaces,
                                     const std::map<localIndex,int>& edgeLocations,
                                     const std::map<localIndex,int>& faceLocations,
                                     const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{
  //do nothing.  We do everything in ProcessNode.
}


//For now, we only enforce the no-lag constraints for the KGD example.
//We need to implement other criteria later.
int Fractunator2D::CheckNodeSplitability( const localIndex nodeID,
                                      NodeManagerT& nodeManager,
                                      FaceManagerT& faceManager,
                                      EdgeManagerT& edgeManager,
                                      const bool prefrac)
{
  //     Return value = -1, this node won't split for sure, don't do any more work;
  //                  = 0, node is a tip, but the fracture connected to it is not saturated yet.  We will only calculate SIF but will not perform splitting.
  //                  = 1, node is a tip and the adjacent fracture is saturated, further check the node;
  //                  = 2, this is a singular node, we need split it.
  //                  = 3, this is an eligible kink, we need to process it as a kink
  int isSplitable = -1;
  const rArray1d* faceFluidPressure = faceManager.GetFieldDataPointer<FieldInfo::pressure>();
  const iArray1d* flowFaceType = faceManager.GetFieldDataPointer<int>("flowFaceType");
  const rArray1d& kinkAngle = nodeManager.GetFieldData<realT>("kinkAngle");

  if (kinkAngle[nodeID] < -0.5 && kinkAngle[nodeID] > -1.5) // This is an interior node with less than 2 external faces.
  {
    isSplitable = -1;
    return (isSplitable);
  }
  else if (kinkAngle[nodeID] < -1.5)// This is a singular node.
  {
    isSplitable = 2;
    return (isSplitable);
  }
  else if (kinkAngle[nodeID] == 0.0) // This is a tip
  {
    if ( faceFluidPressure == NULL || flowFaceType == NULL || m_allowVacuumFrac == 1)  // This is a dry simulation
    {
      isSplitable = 1;
    }
    else
    {

      int nFlowFaceOnNode = 0;
      int nSaturatedFlowFaceOnNode = 0;
      int nExternalFaces = 0;

      // We count the saturated flow cells connected to the tip
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
      {
        if ( (*flowFaceType)[*iface] == 0 )
        {
          nFlowFaceOnNode++;
          if ( (*faceFluidPressure)[*iface] > 0) nSaturatedFlowFaceOnNode++;
        }
      }

      // The number of external faces is counted for the tip.
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
      {
        if ( faceManager.m_isExternal[*iface] == 1) nExternalFaces++;
      }

      if (nExternalFaces >= 4)
      {
        isSplitable = 2;
        // In this case, this node needs to be split regardless of rupture states of faces connected to it.
        return (isSplitable);
      }
      else if  (nFlowFaceOnNode >= 1  )
      {
        isSplitable = 0;
        if (nSaturatedFlowFaceOnNode >=1 ) isSplitable = 1;
      }
    }
  }
  else if (kinkAngle[nodeID]>0 && kinkAngle[nodeID] <= m_maxKinkAngle)
  {
    if ( faceFluidPressure == NULL || flowFaceType == NULL)  // This is a dry simulation
    {
      isSplitable = 3;
    }
    else
    {

      localIndex ancestorNode = nodeID;
      while (nodeManager.m_parentIndex[ancestorNode] < nodeManager.DataLengths())
      {
        ancestorNode = nodeManager.m_parentIndex[ancestorNode];
      } // We nned to do this because when we are checking a child node, we need to see how flow faces are connected to its parent node.


      int nFlowFaceOnNode = 0;
      int nSaturatedFlowFaceOnNode = 0;
      int nExternalFaces = 0;

      // We count the saturated flow cells connected to the ancestor
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[ancestorNode].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[ancestorNode].end() ; ++iface )
      {
        if ( (*flowFaceType)[*iface] == 0 )
        {
          nFlowFaceOnNode++;
          if ( (*faceFluidPressure)[*iface] > 0) nSaturatedFlowFaceOnNode++;
        }
      }

      // The number of external faces is counted for the node.
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
      {
        if ( faceManager.m_isExternal[*iface] == 1) nExternalFaces++;
      }

      if (nExternalFaces != 2 )
      {
        int rank ;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        char msg[200];
        sprintf(msg, "Error! Kink at node %d in partition %d has less than or more than 2 external faces.", int(nodeID), rank);
        throw GPException(msg);

      }

      if (nSaturatedFlowFaceOnNode >=1 ) isSplitable = 3;
    }

  }

  else if (kinkAngle[nodeID]>m_maxKinkAngle && kinkAngle[nodeID] < 181.0)
  {

    //We check for intersection points
    // If there are at least three fracture wings and one of them is saturated, we treat it as a kink although the kink angle is greater than the max kink angle.

    if ( faceFluidPressure == NULL || flowFaceType == NULL)  // This is a dry simulation
    {
      isSplitable = 3;
    }
    else
    {

      localIndex ancestorNode = nodeID;
      while (nodeManager.m_parentIndex[ancestorNode] < nodeManager.DataLengths())
      {
        ancestorNode = nodeManager.m_parentIndex[ancestorNode];
      } // We nned to do this because when we are checking a child node, we need to see how flow faces are connected to its parent node.

      int nFlowFaceOnNode = 0;
      int nSaturatedFlowFaceOnNode = 0;

      // We count the saturated flow cells connected to the ancestor
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[ancestorNode].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[ancestorNode].end() ; ++iface )
      {
        if ( (*flowFaceType)[*iface] == 0 )
        {
          nFlowFaceOnNode++;
          if ( (*faceFluidPressure)[*iface] > 0) nSaturatedFlowFaceOnNode++;
        }
      }


      if (nFlowFaceOnNode>=3 && nSaturatedFlowFaceOnNode>=1) isSplitable=3;
    }

  }



  // If we only care about tips, we don't need to find the parent node
  //    localIndex ancestorNode = nodeID;
  //    while (nodeManager.m_parentIndex[ancestorNode] < nodeManager.DataLengths())
  //    {
  //      ancestorNode = nodeManager.m_parentIndex[ancestorNode];
  //    } // We nned to do this because when we are checking a child node, we need to see how flow faces are connected to its parent node.



  return (isSplitable);

}


int Fractunator2D::CheckEdgeSplitability( const localIndex edgeID,
                                          NodeManagerT& nodeManager,
                                          FaceManagerT& faceManager,
                                          EdgeManagerT& edgeManager,
                                          const bool prefrac)
{
  return(1);
}


realT Fractunator2D::CalculateKinkAngle (const localIndex nodeID,
                                         const NodeManagerT& nodeManager,
                                         const FaceManagerT& faceManager)
{
  localIndex nExternalFaces = 0;
  localIndex faceID[2];
  R1Tensor faceVectorT[2], faceVectorN[2], vec0, vec1;
  realT kinkAngle;

  for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
      iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
  {
    if (faceManager.m_isExternal[*iface] == 1)
    {
      if (nExternalFaces <=1 ) faceID[nExternalFaces] = *iface;
      nExternalFaces++;
    }
  }

  if (nExternalFaces < 2)
  {
    return(-1.0);// throw GPException("Kink angle cannot be calculated for this node. ");
  }
  else if (nExternalFaces >2)
  {
    return (-2.0);
  }

  if (faceManager.m_toNodesRelation[faceID[1]][0] != nodeID)
  {
    // We need to make sure the head of face 0 is the tail of face 1.
    localIndex i = faceID[1];
    faceID[1] = faceID [0];
    faceID[0] = i;
  }

  if (faceManager.m_parentIndex[faceID[0]] < faceManager.DataLengths() && faceManager.m_parentIndex[faceID[0]] == faceManager.m_parentIndex[faceID[1]])
  {
    kinkAngle = 0.0;
  }
  else
  {

    faceVectorT[0] =  (*nodeManager.m_refposition)[faceManager.m_toNodesRelation[faceID[0]][1]];
    faceVectorT[0] -= (*nodeManager.m_refposition)[faceManager.m_toNodesRelation[faceID[0]][0]];
    faceVectorT[0].Normalize();

    faceVectorT[1] =  (*nodeManager.m_refposition)[faceManager.m_toNodesRelation[faceID[1]][1]];
    faceVectorT[1] -= (*nodeManager.m_refposition)[faceManager.m_toNodesRelation[faceID[1]][0]];
    faceVectorT[1].Normalize();

    faceVectorN[0] = faceManager.FaceNormal(nodeManager, faceID[0]);
    faceVectorN[1] = faceManager.FaceNormal(nodeManager, faceID[1]);

    vec0 = faceVectorT[1];
    vec0 -= faceVectorT[0];
    if (vec0.L2_Norm() < 0.0001)
    {
      kinkAngle = 180.0;
    }
    else
    {
      kinkAngle = acos ( 0 - faceVectorT[0] * faceVectorT[1]) / 3.141592653589793238462 * 180.0;

      vec0.Normalize();
      vec1 = faceVectorN[0];
      vec1 += faceVectorN[1];
      vec1 /= 2.0;

      if (vec0 * vec1 < 0) kinkAngle = 360 - kinkAngle;


    }
  }

  return (kinkAngle);

}


int Fractunator2D::ProcessTip (const localIndex nodeID,
                               NodeManagerT& nodeManager,
                               FaceManagerT& faceManager,
                               localIndex splitMode,
                               ModifiedObjectLists& modifiedObjects)
{
  int ival = 0;
  R1Tensor vecTipNorm, vecTip;
  realT SIF = CalculateNodeSIF( nodeID, nodeManager, faceManager, vecTipNorm, vecTip);

  if ( SIF > MinimumToughnessOnNode(nodeID, nodeManager, faceManager) * m_thresholdForEvaluateFace )
  {
    ival = MarkRuptureFaceFromNode ( nodeID, nodeManager, faceManager, vecTipNorm, vecTip, modifiedObjects, splitMode);
  }

  return ival;
}


realT Fractunator2D::CalculateNodeSIF( const localIndex nodeID,
                        NodeManagerT& nodeManager,
                        FaceManagerT& faceManager,
                        R1Tensor& vecTipNorm, R1Tensor& vecTip)
{
  realT rval;
  localIndex nExternalFaces = 0;
  localIndex faceID[2];
  rArray1d& SIF_I = nodeManager.GetFieldData<realT>("SIF_I");
  rArray1d& SIF_II = nodeManager.GetFieldData<realT>("SIF_II");

  SIF_I[nodeID] = 0.0;
  SIF_II[nodeID] = 0.0;


  for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
      iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
  {
    if (faceManager.m_isExternal[*iface] == 1)
    {
      if (nExternalFaces <=1 ) faceID[nExternalFaces] = *iface;
      nExternalFaces++;
    }
  }

  if (nExternalFaces > 2)
  {
    SIF_I[nodeID] = 0.0;
    SIF_II[nodeID] = 0.0;
    throw GPException("Error! This is a singular point, not a tip.  This should not happen!");
  }

  if (faceManager.m_parentIndex[faceID[0]] >= faceManager.DataLengths() || faceManager.m_parentIndex[faceID[0]] != faceManager.m_parentIndex[faceID[1]])
  {
    // Something is wrong! This node has two external faces but they don't share the parent face.
    char msg[200];
    sprintf(msg, "node %d has two external faces but they don't share the parent face.", int(nodeID));
    throw GPException(msg);
  }



  //                           ___________            ^ vecTipNorm
  //                                |                 |
  //                           _____|A'____           |
  //                           _____ ______>O         ----> vecTip
  //                                |A
  //                           _____|_____

  //  Point O is the tip and points A and A' are the two nodes leading to this tip
  //  To make VCCT work for the case where traction is applied along the fracture surface, average nodal force at A and B must be subtracted from that at O.
  //  If A or B happens to be at the boundary of a partition, then the nodal force is not accurate.
  //  Therefore, if we see all elements connected to say, node A are ghosts, we multiply the nodal force calculated by two to correct for this.

  // Which face is A and which is A' used to be arbitrary.  In the new implementation, vecTip X vecTipNorm needs to point out of the plane.

  vecTip = (*nodeManager.m_refposition)[nodeID];
  localIndex j, nodeA[2], faceA, faceAp;

  faceA = faceID[0];
  faceAp = faceID[1];

  j = faceManager.m_parentIndex[faceID[0]];
  j = faceManager.m_toNodesRelation[j][0] + faceManager.m_toNodesRelation[j][1] - nodeID;
  vecTip -= (*nodeManager.m_refposition)[j];
  vecTip.Normalize();

  R1Tensor  cross;
  vecTipNorm = faceManager.FaceNormal(nodeManager, faceA);
  vecTipNorm -= faceManager.FaceNormal(nodeManager, faceAp);
  vecTipNorm.Normalize();

  cross.Cross(vecTip, vecTipNorm);
  if (cross[2] < 0)
  {
    faceA = faceID[1];
    faceAp = faceID[0];
    vecTipNorm *= -1;
  }

  nodeA[0] = faceManager.m_toNodesRelation[faceA][0] + faceManager.m_toNodesRelation[faceA][1] - nodeID;
  nodeA[1] = faceManager.m_toNodesRelation[faceAp][0] + faceManager.m_toNodesRelation[faceAp][1] - nodeID;

  // sanity check
  {
    localIndex ancestor[2];
    for (localIndex i=0; i<2; ++i)
    {
      ancestor[i] = nodeManager.GetParentIndex(nodeA[i]);
    }

    if ( ancestor[0] != ancestor[1])
    {
      std:: cout << nodeID << " " << nodeA[0] << " " << nodeA[1] << " " << ancestor[0] << " " << ancestor[1];
      throw GPException("Error! The two nodes leading to a tip don't share the same parent node." );
    }
  }


  //Calculate nodal force at tip O

  R1Tensor fNodeO, fNodeA[2], distAA, tipOpening;
  localIndex nElemEachSide[2], nGhostElem;


  fNodeO = 0.0;
  nElemEachSide[0] = 0;
  nElemEachSide[1] = 0;

  for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
      k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
  {
    ElementRegionT* elementRegion = k->first;
    localIndex iEle = k->second;
    R1Tensor fN;

    elementRegion->CalculateNodalForcesFromOneElement( nodeID, iEle, nodeManager, fN);

    R1Tensor xEle;
    const localIndex* elemToNodeMap = elementRegion->m_toNodesRelation[iEle];

    for (localIndex i = 0; i < (elementRegion->m_numNodesPerElem); ++i)
    {
      xEle += (*nodeManager.m_refposition)[elemToNodeMap[i]];
    }
    xEle /= elementRegion->m_numNodesPerElem;
    xEle -= (*nodeManager.m_refposition)[nodeID];

    if (Dot(xEle, vecTipNorm) > 0)
    {
      nElemEachSide[0] += 1;
      fNodeO += fN;
    }
    else
    {
      nElemEachSide[1] +=1;
      fNodeO -= fN;
    }
  }

  if (nElemEachSide[0]>=1 && nElemEachSide[1]>=1) fNodeO /= 2.0; //We have contributions from both sides.

  // We have to subtract the nodal force at nodes A and A' to take into account the effects of surface traction along the fracture.
  for (localIndex i=0; i<2; ++i)
  {
    nGhostElem = 0;
    fNodeA[i] = 0.0;
    for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeA[i]].begin() ;
        k!=nodeManager.m_toElementsRelation[nodeA[i]].end() ; ++k )
    {
      ElementRegionT* elementRegion = k->first;
      localIndex iEle = k->second;
      iArray1d& elem_is_ghost = elementRegion->GetFieldData<FieldInfo::ghostRank>();
      R1Tensor fN;

      elementRegion->CalculateNodalForcesFromOneElement( nodeA[i], iEle, nodeManager, fN);

      fNodeA[i] += fN;

      if (elem_is_ghost[iEle] >= 0) nGhostElem +=1;

    }
    if (nGhostElem == nodeManager.m_toElementsRelation.size())
    {
      // All elements connected to this node are ghost elements.  This implies that half of the elements are in the next partition.
      fNodeA[i] *= 2.0;
    }

  }

  realT tipForce[2];

  tipForce[0] = Dot(fNodeO, vecTipNorm) + Dot(fNodeA[0], vecTipNorm) / 2.0 - Dot(fNodeA[1], vecTipNorm) /2.0;
  tipForce[1] = Dot(fNodeO, vecTip) + Dot(fNodeA[0], vecTip) / 2.0 - Dot(fNodeA[1], vecTip) /2.0;

  distAA = (*nodeManager.m_displacement)[nodeA[1]];
  distAA -= (*nodeManager.m_displacement)[nodeA[0]];

  tipOpening[0] = Dot(distAA, vecTipNorm);
  tipOpening[1] = Dot(distAA, vecTip);

  realT tipArea;
  tipArea = Dot((*nodeManager.m_refposition)[nodeA[0]] - (*nodeManager.m_refposition)[nodeID],
                (*nodeManager.m_refposition)[nodeA[0]] - (*nodeManager.m_refposition)[nodeID]);
  tipArea = pow(tipArea, 0.5);

  SIF_I[nodeID] = pow(fabs(tipForce[0] * tipOpening[0] / 2.0 / tipArea), 0.5);
  SIF_II[nodeID] = pow(fabs(tipForce[1] * tipOpening[1] / 2.0 / tipArea), 0.5);

  if (tipForce[1] < 0) SIF_II[nodeID] *= -1;
  if (tipOpening[0] < 0) SIF_I[nodeID] *= -1;

  if (SIF_I[nodeID] > 0.0)
  {
    rval = pow(SIF_I[nodeID]*SIF_I[nodeID]+SIF_II[nodeID]*SIF_II[nodeID], 0.5);
  }
  else
  {
    rval = -1.0;
  }

  return rval;
}


int Fractunator2D::MarkRuptureFaceFromNode ( const localIndex nodeID,
                                              NodeManagerT& nodeManager,
                                              FaceManagerT& faceManager,
                                              R1Tensor& vecTipNorm,
                                              R1Tensor& vecTip,
                                              ModifiedObjectLists& modifiedObjects,
                                              const int splitMode)
{

  // We only deal with SIF-based criterion here.  The stress-based criterion is handled in the UpdateRuptureState function
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  rArray1d& SIFonFace = faceManager.GetFieldData<realT>("SIFonFace");
  iArray1d& isFaceSeparable = faceManager.GetFieldData<int>("isSeparable");
  rArray1d& faceToughness = faceManager.GetFieldData<realT>("faceToughness");
  rArray1d& SIF_I = nodeManager.GetFieldData<realT>("SIF_I");
  rArray1d& SIF_II = nodeManager.GetFieldData<realT>("SIF_II");

  localIndex maxSIFFace, nEffFace(0);
  realT maxSIF = std::numeric_limits<realT>::min();

  for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
      iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
  {
    R1Tensor vecFace;
    localIndex node1 = faceManager.m_toNodesRelation[*iface][0] + faceManager.m_toNodesRelation[*iface][1] - nodeID;

    vecFace = (*nodeManager.m_refposition)[node1];
    vecFace -= (*nodeManager.m_refposition)[nodeID];
    vecFace.Normalize();

    if (Dot(vecTip, vecFace) > cos(m_maxTurnAngle))
    {
      realT thetaFace = acos(Dot(vecTip,vecFace)*0.999999);  // We multiply this by 0.9999999 to avoid an exception caused by acos a number slightly larger than 1.

      if (Cross(vecTip,vecFace)[2] < 0.0)
      {
        thetaFace *= -1.0;
      }

      SIFonFace[*iface] = cos(thetaFace / 2.0) * ( SIF_I[nodeID] * cos(thetaFace / 2.0) * cos(thetaFace / 2.0) - 1.5 * SIF_II[nodeID] * sin(thetaFace) );

      // Hopefully we won't have to do the mixed voting thing as in 3D.
      realT SIFRatio = SIFonFace[*iface] / faceToughness[*iface];
      if (SIFRatio > 1.0 && isFaceSeparable[*iface] == 1)
      {
        nEffFace++;
        if (SIFRatio > maxSIF)
        {
          maxSIF = SIFRatio;
          maxSIFFace = *iface;
        }
      }
    }
  }

  if (nEffFace >= 1 )
  {
    ruptureState[maxSIFFace] = 1;
    modifiedObjects.modifiedFaces.insert(maxSIFFace);

    return 1;
  }
  else
  {
    return 0;
  }
}



int Fractunator2D::ProcessKink (const localIndex nodeID,
                                NodeManagerT& nodeManager,
                                FaceManagerT& faceManager,
                                ModifiedObjectLists& modifiedObjects)
{

  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  rArray1d& stressNOnFace = faceManager.GetFieldData<realT>("stressNOnFace");
  iArray1d& isFaceSeparable = faceManager.GetFieldData<int>("isSeparable");

  localIndex maxTensionFace, nEffFace;
  realT maxTension = std::numeric_limits<realT>::min();
  nEffFace = 0;
  for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
      iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface )
  {

    if (faceManager.m_toElementsRelation[*iface].size() == 2  &&
        stressNOnFace[*iface] > m_kinkStrength &&
        CheckOrphanElement(faceManager, *iface) == 0 &&
        isFaceSeparable[*iface] == 1)
    {
      nEffFace++;
//      ruptureState[*iface] = 0;
      if (nEffFace == 1)
      {
        maxTensionFace = *iface;
        maxTension = stressNOnFace[maxTensionFace];
      }
      else if ( stressNOnFace[*iface] > maxTension)
      {
        maxTensionFace = *iface;
        maxTension =stressNOnFace[maxTensionFace];
      }

    }
  }

  if (nEffFace >= 1)
  {
    ruptureState[maxTensionFace] = 1;
    modifiedObjects.modifiedFaces.insert(maxTensionFace);
    return 1;
  }
  else
  {
    return 0;
  }

}

void Fractunator2D::WriteTopology(ofstream& outfile_,
                                  NodeManagerT& nodeManager, ElementManagerT& elementManager,
                                  FaceManagerT& faceManager, EdgeManagerT& edgeManager)
{



  outfile_ << "Node to Element map" << std::endl;
  for (localIndex i = 0; i != nodeManager.DataLengths(); ++i)
  {
    outfile_ << "node " << setw(5) << i << ": ";
    for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[i].begin() ;
        k!=nodeManager.m_toElementsRelation[i].end() ; ++k )
    {
      localIndex iEle = k->second;
      outfile_ << setw(8) << iEle << ", ";
    }
    outfile_ << std::endl;
  }


  outfile_ << "Node to Face map" << std::endl;
  for (localIndex i = 0; i != nodeManager.DataLengths(); ++i)
  {
    outfile_ << "node " << setw(5) << i << ": ";

    for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[i].begin() ;
        iface!=nodeManager.m_nodeToFaceMap[i].end() ; ++iface )
    {
      outfile_ << setw(8) << *iface << ", ";
    }
    outfile_ << std::endl;
  }


  outfile_ << "Node to Edge map" << std::endl;
  for (localIndex i = 0; i != nodeManager.DataLengths(); ++i)
  {
    outfile_ << "node " << setw(5) << i << ": ";

    for( lSet::const_iterator iedge=nodeManager.m_nodeToEdgeMap[i].begin() ;
        iedge!=nodeManager.m_nodeToEdgeMap[i].end() ; ++iedge )
    {
      outfile_ << setw(8) << *iedge << ", ";
    }
    outfile_ << std::endl;
  }


  outfile_ << "Element to Node map" << std::endl;
  const ElementRegionT& elemRegion = (elementManager.m_ElementRegions.begin())->second;

  for (localIndex i = 0; i != elemRegion.DataLengths(); ++i)
  {
    outfile_ << "element " << setw(5) << i << ": ";

    for (localIndex j = 0; j < elemRegion.m_toNodesRelation.Dimension(1); ++j)
    {
      outfile_ << setw(8) << elemRegion.m_toNodesRelation[i][j] << ", ";
    }
    outfile_ << std::endl;
  }


  outfile_ << "Element to Face map" << std::endl;

  for (localIndex i = 0; i != elemRegion.DataLengths(); ++i)
  {
    outfile_ << "element " << setw(5) << i << ": ";

    for (localIndex j = 0; j < elemRegion.m_toFacesRelation.Dimension(1); ++j)
    {
      outfile_ << setw(8) << elemRegion.m_toFacesRelation[i][j] << ", ";
    }
    outfile_ << std::endl;
  }


  //  outfile_ << "Element to Edge map" << std::endl;
  //
  //  for (localIndex i = 0; i != elemRegion.DataLengths(); ++i)
  //  {
  //    outfile_ << "element " << setw(5) << i << ": ";
  //
  //    for (localIndex j = 0; j < elemRegion.m_toEdgesRelation.Dimension(1); ++j)
  //    {
  //      outfile_ << setw(8) << elemRegion.m_toEdgesRelation[i][j] << ", ";
  //    }
  //    outfile_ << std::endl;
  //  }


  outfile_ << "Face to Node map" << std::endl;

  for (localIndex i = 0; i != faceManager.DataLengths(); ++i)
  {
    outfile_ << "face " << setw(5) << i << ": ";

    for( lArray1d::iterator j = faceManager.m_toNodesRelation[i].begin() ;
        j!=faceManager.m_toNodesRelation[i].end() ; ++j )
    {
      outfile_ << setw(8) << *j << ", ";
    }
    outfile_ << std::endl;
  }

  outfile_ << "Face to Edge map" << std::endl;

  for (localIndex i = 0; i != faceManager.DataLengths(); ++i)
  {
    outfile_ << "face " << setw(5) << i << ": ";

    for( lArray1d::iterator j = faceManager.m_toEdgesRelation[i].begin() ;
        j!=faceManager.m_toEdgesRelation[i].end() ; ++j )
    {
      outfile_ << setw(8) << *j << ", ";
    }
    outfile_ << std::endl;
  }

  outfile_ << "Face to Element map" << std::endl;

  for (localIndex i = 0; i != faceManager.DataLengths(); ++i)
  {
    outfile_ << "face " << setw(5) << i << ": ";

    for (Array1dT< std::pair< ElementRegionT*, localIndex > >::iterator iter = faceManager.m_toElementsRelation[i].begin();
        iter != faceManager.m_toElementsRelation[i].end(); ++iter )
    {
      //      ElementRegionT* elementRegion = iter->first;
      localIndex iEle= iter->second;
      outfile_ << setw(8) << iEle << ", ";
    }
    outfile_ << std::endl;
  }


  outfile_ << "Edge to Node map" << std::endl;

  for (localIndex i = 0; i != edgeManager.DataLengths(); ++i)
  {
    outfile_ << "edge " << setw(5) << i << ": ";

    for( unsigned int a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a )
    {
      outfile_ << setw(8) << edgeManager.m_toNodesRelation(i,a) << ", ";
    }
    outfile_ << std::endl;
  }


  outfile_ << "Edge to Face map" << std::endl;

  for (localIndex i = 0; i != edgeManager.DataLengths(); ++i)
  {
    outfile_ << "edge " << setw(5) << i << ": ";

    for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[i].begin() ;
        iface!=edgeManager.m_toFacesRelation[i].end() ; ++iface )
    {
      outfile_ << setw(8) << *iface << ", ";
    }
    outfile_ << std::endl;
  }




  outfile_.close();


}

/// Register solver in the solver factory
REGISTER_FRACTUNATOR( Fractunator2D)


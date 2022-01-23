/**
 * @file Fractunator.h
 * @author settgast1
 * @date Jul 14, 2011
 */

#ifndef FRACTUNATOR2_H_
#define FRACTUNATOR2_H_


#include "SurfaceGeneration/FractunatorBase.h"

class Fractunator2 : public FractunatorBase
{
public:
  Fractunator2();
  virtual ~Fractunator2();

  void PreexistingFracture2D( NodeManagerT& nodeManager,
                   EdgeManagerT& edgeManager,
                   FaceManagerT& faceManager,
                   ExternalFaceManagerT& externalFaceManager,
                   ElementManagerT& elementManager,
                   SpatialPartition& partition,
                   const bool prefrac);

  static std::string FractunatorName() { return "Fractunator2"; }

private:



  bool FindFracturePlanes( const localIndex nodeID,
                           const NodeManagerT& nodeManager,
                           const EdgeManagerT& edgeManager,
                           const FaceManagerT& faceManager,
                           const Array1dT<lSet>& nodesToRupturedFaces,
                           const Array1dT<lSet>& edgesToRupturedFaces,
                           lSet& separationPathFaces,
                           std::map<localIndex,int>& edgeLocations,
                           std::map<localIndex,int>& faceLocations,
                           std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations  );



  void PerformFracture( const localIndex nodeID,
                        NodeManagerT& nodeManager,
                        EdgeManagerT& edgeManager,
                        FaceManagerT& faceManager,
                        ExternalFaceManagerT& externalFaceManager,
                        ElementManagerT& elementManager,
                        ModifiedObjectLists& modifiedObjects,
                        Array1dT<lSet>& nodesToRupturedFaces,
                        Array1dT<lSet>& edgesToRupturedFaces,
                        const lSet& separationPathFaces,
                        const std::map<localIndex,int>& edgeLocations,
                        const std::map<localIndex,int>& faceLocations,
                        const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations );

bool SetLocations( const int side,
                   const lSet& separationPathFaces,
                   const FaceManagerT& faceManager,
                   const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                   std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                   std::map<localIndex,int>& edgeLocations,
                   std::map<localIndex,int>& faceLocations,
                   std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations );


};



#endif /* FRACTUNATOR_H_ */

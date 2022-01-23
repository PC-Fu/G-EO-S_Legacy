/**
 * @file Fractunator.h
 * @author settgast1
 * @date Jul 14, 2011
 */

#ifndef FRACTUNATOR3_H_
#define FRACTUNATOR3_H_

#include "SurfaceGeneration/FractunatorBase.h"


class SpatialPartition;
class ModifiedObjectLists;

class Fractunator3 : public FractunatorBase
{
public:
  Fractunator3();
  virtual ~Fractunator3();

  void Initialize( NodeManagerT& nodeManager,
                   EdgeManagerT& edgeManager,
                   FaceManagerT& faceManager,
                   ElementManagerT& elementManager);

  static std::string FractunatorName() { return "Fractunator3"; }

  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  void RegisterFieldsAndMaps( NodeManagerT& nodeManager,
                              EdgeManagerT& edgeManager,
                              FaceManagerT& faceManager );

  void PreexistingFracture2D( NodeManagerT& nodeManager,
                   EdgeManagerT& edgeManager,
                   FaceManagerT& faceManager,
                   ExternalFaceManagerT& externalFaceManager,
                   ElementManagerT& elementManager,
                   SpatialPartition& partition,
                   const bool prefrac);

  int SeparationDriver( NodeManagerT& nodeManager,
                         EdgeManagerT& edgeManager,
                         FaceManagerT& faceManager,
                         ExternalFaceManagerT& externalFaceManager,
                         ElementManagerT& elementManager,
                         SpatialPartition& partition,
                         const bool prefrac,
                         const realT time);

  void IdentifyRupturedFaces( NodeManagerT& nodeManager,
                              EdgeManagerT& edgeManager,
                              FaceManagerT& faceManager,
                              ElementManagerT& elementManager,
                              SpatialPartition& partition,
                              const bool prefrac  );

  realT CalculateEdgeSIF ( const localIndex edgeID,
                         NodeManagerT& nodeManager,
                         EdgeManagerT& edgeManager,
                         FaceManagerT& faceManager,
                         ElementManagerT& elementManager,
                         R1Tensor& vecTipNorm,
                         R1Tensor& vecTip);

  void MarkRuptureFaceFromEdge ( const localIndex edgeID,
                                 NodeManagerT& nodeManager,
                                 EdgeManagerT& edgeManager,
                                 FaceManagerT& faceManager,
                                 R1Tensor& vecTipNorm,
                                 R1Tensor& vecTip,
                                 ModifiedObjectLists& modifiedObjects,
                                 const int edgeMode);



//  FaceManagerT m_virtualFaces;

  EdgeManagerT m_virtualEdges;

  NodeManagerT m_virtualNodes;

  rArray1d m_insituStress3D;  // components 11, 22, 33, 23, 31, 12


  int m_markExtendedLayer;

  realT m_saturationPressureCuttoff;

  realT m_faceToEdgeProjectionTol;

private:
  void UpdateRuptureStates( NodeManagerT& nodeManager,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager,
                            ElementManagerT& elementManager,
                            Array1dT<lSet>& nodesToRupturedFaces,
                            Array1dT<lSet>& edgesToRupturedFaces,
                            const bool prefrac = false );
  void PostUpdateRuptureStates( NodeManagerT& nodeManager,
                                EdgeManagerT& edgeManager,
                                FaceManagerT& faceManager,
                                ElementManagerT& elementManager,
                                Array1dT<lSet>& nodesToRupturedFaces,
                                Array1dT<lSet>& edgesToRupturedFaces);

  int CheckEdgeSplitability( const localIndex edgeID,
                                NodeManagerT& nodeManager,
                                FaceManagerT& faceManager,
                                EdgeManagerT& edgeManager,
                                const bool prefrac);

  int CheckNodeSplitability( const localIndex nodeID,
                                NodeManagerT& nodeManager,
                                FaceManagerT& faceManager,
                                EdgeManagerT& edgeManager,
                                const bool prefrac);

  void UpdatePathCheckingArrays();

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


bool SetLocations( const lSet& separationPathFaces,
                   const FaceManagerT& faceManager,
                   const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                   const std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                   std::map<localIndex,int>& edgeLocations,
                   std::map<localIndex,int>& faceLocations,
                   std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations );

bool SetElemLocations( const int side,
                       const std::pair< ElementRegionT*, localIndex >& elem,
                       const lSet& separationPathFaces,
                       const FaceManagerT& faceManager,
                       const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                       const std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                       std::map<localIndex,int>& edgeLocations,
                       std::map<localIndex,int>& faceLocations,
                       std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations );



  void ApplyGapDamping( NodeManagerT& nodeManager,
                        const FaceManagerT& faceManager,
                        const realT dt  );

  virtual void WriteSiloDerived( SiloFile& siloFile,
                                 const int cycleNum,
                                 const realT problemTime,
                                 const bool isRestart );

  virtual void ReadSiloDerived( const SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart );


  realT CalculateKinkAngle (const localIndex nodeID,
                            const NodeManagerT& nodeManager,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager);

  void CalculateKinkAngles (FaceManagerT& faceManager,
                            EdgeManagerT& edgeManager,
                            NodeManagerT& nodeManager,
                            ModifiedObjectLists& modifiedObjects,
                            const bool prefrac);

};



#endif /* FRACTUNATOR3_H_ */

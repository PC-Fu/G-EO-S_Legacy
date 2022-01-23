/**
 * @file Fractunator2D.h
 * @author fu4
 * @date Oct. 29, 2012
 */

#ifndef FRACTUNATOR2D_H_
#define FRACTUNATOR2D_H_

#include "SurfaceGeneration/FractunatorBase.h"
#include "Common/Common.h"
//#include "ObjectManagers/PhysicalDomainT.h"
//#include "IO/ticpp/TinyXMLParser.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;


class SpatialPartition;
//class ModifiedObjectLists;

class Fractunator2D : public FractunatorBase
{
public:
  Fractunator2D();
  virtual ~Fractunator2D();

  void Initialize( NodeManagerT& nodeManager,
                   EdgeManagerT& edgeManager,
                   FaceManagerT& faceManager,
                   ElementManagerT& elementManager);


  static std::string FractunatorName() { return "Fractunator2D"; }


  void RegisterFieldsAndMaps( NodeManagerT& nodeManager,
                              EdgeManagerT& edgeManager,
                              FaceManagerT& faceManager );

  virtual int SeparationDriver( NodeManagerT& nodeManager,
                                EdgeManagerT& edgeManager,
                                FaceManagerT& faceManager,
                                ExternalFaceManagerT& externalFaceManager,
                                ElementManagerT& elementManager,
                                SpatialPartition& partition,
                                const bool prefrac,
                                const realT time);

  void ReadXML( TICPP::HierarchicalDataNode& hdn );



  void WriteTopology(ofstream& outfile_,
                     NodeManagerT& nodeManager, ElementManagerT& elementManager,
                     FaceManagerT& faceManager, EdgeManagerT& edgeManager);

  void PreexistingFracture2D( NodeManagerT& nodeManager,
                   EdgeManagerT& edgeManager,
                   FaceManagerT& faceManager,
                   ExternalFaceManagerT& externalFaceManager,
                   ElementManagerT& elementManager,
                   SpatialPartition& partition,
                   const bool prefrac);



  FaceManagerT m_virtualFaces;
  //  void WriteSilo( SiloFile& restartFile ) const;
  //
  //  void ReadSilo( const SiloFile& restartFile );



private:


  void TrackAlongAFractures (FaceManagerT& faceManager,
                             NodeManagerT& nodeManager,
                             R1Tensor x0,
                             R1Tensor x1,
                             ModifiedObjectLists& modifiedObjects);

  void TrackAlongShortestPath (FaceManagerT& faceManager,
                               NodeManagerT& nodeManager,
                               R1Tensor x0,
                               R1Tensor x1,
                               ModifiedObjectLists& modifiedObjects);


  void CreatePreexistingFractures (NodeManagerT& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager,
                                   ExternalFaceManagerT& externalFaceManager,
                                   ElementManagerT& elementManager,
                                   SpatialPartition& partition,
                                   const bool prefrac);

  void IdentifyRupturedFaces( NodeManagerT& nodeManager,
                              EdgeManagerT& edgeManager,
                              FaceManagerT& faceManager,
                              ElementManagerT& elementManager,
                              SpatialPartition& partition,
                              const bool prefrac  );

  void FindClosestNode (NodeManagerT& nodeManager,
                        R1Tensor& xPt,
                        localIndex& iNdMinDist );

  realT CalculateKinkAngle (const localIndex nodeID,
                            const NodeManagerT& nodeManager,
                            const FaceManagerT& faceManager);

  bool EvaluateAndSplitNode( const localIndex nodeID,
                            NodeManagerT& nodeManager,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager,
                            ExternalFaceManagerT& externalFaceManager,
                            ElementManagerT& elementManager,
                            ModifiedObjectLists& modifiedObjects);

//  bool CheckAndSplitNode( const localIndex nodeID,
//                            NodeManagerT& nodeManager,
//                            EdgeManagerT& edgeManager,
//                            FaceManagerT& faceManager,
//                            ExternalFaceManagerT& externalFaceManager,
//                            ElementManagerT& elementManager,
//                            ModifiedObjectLists& modifiedObjects,
//                            const bool prefrac);

  int CheckNodeSplitability( const localIndex nodeID,
                             NodeManagerT& nodeManager,
                             FaceManagerT& faceManager,
                             EdgeManagerT& edgeManager,
                             const bool prefrac);

  int CheckEdgeSplitability( const localIndex edgeID,
                             NodeManagerT& nodeManager,
                             FaceManagerT& faceManager,
                             EdgeManagerT& edgeManager,
                             const bool prefrac);

  int MakeDividingVectors (localIndex nodeID,
                           R1Tensor& vecWingNorm,
                           R1Tensor vecWings[],
                           localIndex dividingFaces[],
                           NodeManagerT& nodeManager,
                           FaceManagerT& faceManager);

  int MySideAlongFracture (localIndex nodeID,
                           const std::pair< ElementRegionT*, localIndex >& elem,
                           int iDivType,
                           R1Tensor vecWingNorm,
                           R1Tensor vecWings[],
                           NodeManagerT& nodeManager,
                           ElementManagerT& elementManager );

  int ProcessTip (const localIndex nodeID,
                  NodeManagerT& nodeManager,
                  FaceManagerT& faceManager,
                  localIndex splitMode,
                  ModifiedObjectLists& modifiedObjects);

  realT CalculateNodeSIF( const localIndex nodeID,
                          NodeManagerT& nodeManager,
                          FaceManagerT& faceManager,
                          R1Tensor& vecTipNorm, R1Tensor& vecTip);

  int MarkRuptureFaceFromNode ( const localIndex nodeID,
                                 NodeManagerT& nodeManager,
                                 FaceManagerT& faceManager,
                                 R1Tensor& vecTipNorm,
                                 R1Tensor& vecTip,
                                 ModifiedObjectLists& modifiedObjects,
                                 const int splitMode);


  int ProcessKink (const localIndex nodeID,
                  NodeManagerT& nodeManager,
                  FaceManagerT& faceManager,
                  ModifiedObjectLists& modifiedObjects);

  void MarkPreexistingFractures (FaceManagerT& faceManager,
                                 NodeManagerT& nodeManager,
                                 SpatialPartition& partition,
                                 ModifiedObjectLists& modifiedObjects);


  Array1dT< R1Tensor > m_x0_PreFrac;
  Array1dT< R1Tensor > m_x1_PreFrac;
  R1Tensor m_insituStress_2D;

  realT m_maxKinkAngle;
  realT m_kinkStrength;

  //std::ofstream outfile;


  //std::string m_separableSet;
  //realT m_failgap;

  //  virtual void UpdateCohesiveZones( NodeManagerT& nodeManager,
  //                                   EdgeManagerT& edgeManager,
  //                                   FaceManagerT& faceManager ){}


  //protected:
  //  void UpdateRuptureStates( NodeManagerT& nodeManager,
  //                            EdgeManagerT& edgeManager,
  //                            FaceManagerT& faceManager,
  //                            ElementManagerT& elementManager,
  //                            Array1dT<lSet>& nodesToRupturedFaces,
  //                            Array1dT<lSet>& edgesToRupturedFaces );
  //
  //
  //private:
  //
  //  bool ProcessNode( const localIndex nodeID,
  //                            NodeManagerT& nodeManager,
  //                            EdgeManagerT& edgeManager,
  //                            FaceManagerT& faceManager,
  //                            Array1dT<lSet>& nodesToRupturedFaces,
  //                            Array1dT<lSet>& edgesToRupturedFaces,
  //                            ElementManagerT& elementManager,
  //                            ModifiedObjectLists& modifiedObjects );
  //
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
};



#endif /* FRACTUNATOR2D_H_ */

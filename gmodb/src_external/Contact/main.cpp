#include "CommonPlaneContact.h"

bool TestFaceNodeContributions1()
{
  //make 2 faces perpendicular to the x-axis
  const int nfaces = 2;
  const int nnodes = 3*nfaces;

  //nodes on face 0 are at x = 0 w/ nx = -1; nodes on face 1 are at x = 0 with nx = 1
  const double xNodes[] = {0,0.1,0,0,0,-0.1,0,0,0.1,0,0,0.1,0,0,-0.1,0,0.1,0};
  //nodes on face 0 are at dx = 1; nodes on face 1 are at dx = -1 (increasing overlap)
  const double dxNodes[] = {-1,0,0,-1,0,0,-1,0,0,1,0,0,1,0,0,1,0,0};
  //no initial accelerations
  double ddxNodes[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  //unit nodal masses
  const double mNodes[] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};

  int faceNodeCounts[] = {3,3};
  int faceOffsets[] = {0,3};
  int faceNodes[] = {0,1,2,3,4,5};
  int faceIDs[] = {0,1};
  int faceGhosts[] = {0,0};
  int facePairCounts[] = {1,0};
  int facePairOffsets[] = {0,1};

  const double kbulk[] = {1.e10, 1.e10};
  const double kshear[] = {1.e8, 1.e8};
  const double viscosity[] = {0.01, 0.01};
  const double cohesion[] = {0, 0};
  const double cohesionLimit[] = {1.e11, 1.e11};
  const double frictionSlope[] = {0.6, 0.6};
  const double residualSlope[] = {0.3, 0.3};
  const double criticalShearSlip[] = {1e-4, 1e-4};
  const double dilation[] = {0., 0.};
  double dilationStress[] = {0., 0.};
  const double dilationAngle[] = {0.7, 0.7};
  const double dilationLimit[] = {1.e10, 1.e10};
  const double tensileSoftening[] = {0., 0.};
  const double softening[] = {0., 0.};
  double tensileDamage[] = {0., 0.};
  const double tensileCutOff[] = {1.e10, 1.e10};
  const double pf[] = {0., 0.};

  const int face2[] = {1};
  int active[] = {0};
  double aperture[] = {0.};
  double normalApproach[] = {0};
  double normalApproachMax[] = {0};
  double normal[] = {1,0,0};
  double forceNormal[] = {0};
  double forceDilation[] = {0};
  double forceTangential[] = {0};
  double slip[] = {0};

  //polygons
  //double xpoly[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  std::vector<double> xpoly;
  int npoly;
  int* poly_facepairs = new int[facePairOffsets[nfaces-1] + facePairCounts[nfaces-1]];
  int* poly_offsets = new int[facePairOffsets[nfaces-1] + facePairCounts[nfaces-1]];
  int* poly_counts = new int[facePairOffsets[nfaces-1] + facePairCounts[nfaces-1]];

  const double dt = 1e-6;
  {
    CommonPlaneContact::FaceNodeContributions (
        nnodes,
        xNodes,
        dxNodes,
        ddxNodes,
        mNodes,

        nfaces,
        faceNodeCounts,
        faceOffsets,
        facePairCounts,
        facePairOffsets,
        faceNodes,
        faceIDs,
        faceGhosts,

        kbulk,
        kshear,
        viscosity,
        cohesion,
        cohesionLimit,
        frictionSlope,
        residualSlope,
        criticalShearSlip,
        dilation,
        dilationStress,
        dilationAngle,
        dilationLimit,
        tensileSoftening,
        softening,
        tensileDamage,
        tensileCutOff,
        pf,

        face2,
        active,
        aperture,
        normalApproach,
        normalApproachMax,
        normal,
        forceNormal,
        forceDilation,
        forceTangential,
        slip,

        npoly,
        poly_facepairs,
        poly_offsets,
        poly_counts,
        xpoly,

        dt);
  }

  delete[] poly_facepairs;
  delete[] poly_offsets;
  delete[] poly_counts;
  return true;
}

bool TestPolygonalAreaOfIntersection3_1()
{
  double a1[] = {0,1.1,-0.9,0,0,-1.5,1.5,-0.5,1,0.5};
  unsigned int num1 = 5;
  double a2[] = {-0.9,-1,1.25,-1.1,1,0,0,1,-1,0,-1.25,-0.5};
  unsigned int num2 = 6;

  //double* aa = new double[4*(num1 > num2 ? num1 : num2)];
  std::vector<double> aa;
  //aa.resize(4*(num1 > num2 ? num1 : num2), 0);
  unsigned int npoly = 0;
  double contactArea = CommonPlaneContact::PolygonalAreaOfIntersection3(a1, a2, num1, num2, aa, npoly);
  double contactArea1 = CommonPlaneContact::PolygonalAreaOfIntersection2(a1, a2, num1, num2);
  return fabs(contactArea - contactArea1) < 1e-10 && npoly == 7;
  //delete[] aa;
}

bool TestPolygonalAreaOfIntersection3_2()
{
  double a1[] = {1,0,0,1,-1,0,0,-1};
  unsigned int num1 = 4;
  double a2[] = {1.001,0,0,1.001,-1.001,0,0,-1.001};
  unsigned int num2 = 4;

  //double* aa = new double[4*(num1 > num2 ? num1 : num2)];
  std::vector<double> aa;
  //aa.resize(4*(num1 > num2 ? num1 : num2), 0);
  unsigned int npoly;
  double contactArea = CommonPlaneContact::PolygonalAreaOfIntersection3(a1, a2, num1, num2, aa, npoly);
  double contactArea1 = CommonPlaneContact::PolygonalAreaOfIntersection2(a1, a2, num1, num2);
  return fabs(contactArea - contactArea1) < 1e-10 && npoly == 4;
  //delete[] aa;
}

bool TestPolygonalAreaOfIntersection3_3()
{
  double a1[] = {1,0,0,1.5,-1,0,0,-1};
  unsigned int num1 = 4;
  double a2[] = {1.001,0,0,1.001,-1.001,0,0,-1.001};
  unsigned int num2 = 4;

  //double* aa = new double[4*(num1 > num2 ? num1 : num2)];
  std::vector<double> aa;
  //aa.resize(4*(num1 > num2 ? num1 : num2), 0);
  unsigned int npoly;
  double contactArea = CommonPlaneContact::PolygonalAreaOfIntersection3(a1, a2, num1, num2, aa, npoly);
  double contactArea1 = CommonPlaneContact::PolygonalAreaOfIntersection2(a1, a2, num1, num2);
  return fabs(contactArea - contactArea1) < 1e-10 && npoly == 6;
  //delete[] aa;
}

bool TestPointInTriangle1()
{
  double uv[] = {-1,0,1,0,0,1};

  double uv0[] = {-0.999,1.e-6};
  double uv1[] = {-1.0001,1.e-6};

  double uv2[] = {0.999,1.e-6};
  double uv3[] = {1.0001,1.e-6};

  double uv4[] = {0,0.9999};
  double uv5[] = {0,1.0001};

  int in[] = {-1,-1,-1,-1,-1,-1};

  in[0] = CommonPlaneContact::PointInTriangle(uv, uv0);
  in[1] = CommonPlaneContact::PointInTriangle(uv, uv1);
  in[2] = CommonPlaneContact::PointInTriangle(uv, uv2);
  in[3] = CommonPlaneContact::PointInTriangle(uv, uv3);
  in[4] = CommonPlaneContact::PointInTriangle(uv, uv4);
  in[5] = CommonPlaneContact::PointInTriangle(uv, uv5);

  int in0 = 0;
  for(int i = 0; i < 6; ++i)
    in0 += in[i];
  return in0 == 3;
}

int main(int* argc, char** argv)
{
  bool t0 = TestPointInTriangle1();
  bool t1 = TestPolygonalAreaOfIntersection3_1();
  bool t2 = TestPolygonalAreaOfIntersection3_2();
  bool t3 = TestPolygonalAreaOfIntersection3_3();
  bool t4 = TestFaceNodeContributions1();
  return 0;
}

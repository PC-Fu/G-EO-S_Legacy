

ARRAY_OBJS ?=\
	src/ArrayT/bufvector.o

BC_OBJS ?=\
	src/BoundaryConditions/BoundaryConditions.o \
	src/BoundaryConditions/ApplyBoundaryConditions.o 

CONTACT_OBJS ?=\
	src/Contact/SpatialSorterFactory.o \
	src/Contact/SpatialSorterBase.o \
	src/Contact/CellVerlet.o \
	src/Contact/NSquared.o 


DATASTRUCTURES_OBJS ?=\
	src/DataStructures/VectorFields/ElementRegionT.o \
	src/DataStructures/VectorFields/NodeManagerT.o \
	src/DataStructures/VectorFields/ObjectDataStructureBaseT.o \
	src/DataStructures/VectorFields/TempODS.o \
	src/DataStructures/EncapsulatedFields/ElementRegion.o \
	src/DataStructures/EncapsulatedFields/FieldRegistry.o \
	src/DataStructures/EncapsulatedFields/NodeManager.o \
	src/DataStructures/EncapsulatedFields/ObjectManager.o \
	src/DataStructures/InterObjectRelation.o


ELEMENT_OBJS ?=\
	src/ElementLibrary/FiniteElementBase.o \
	src/ElementLibrary/Polynomial.o \
	src/ElementLibrary/Basis.o \
	src/ElementLibrary/LagrangeBasis.o \
	src/ElementLibrary/Quadrature.o \
	src/ElementLibrary/GaussQuadrature.o \
	src/ElementLibrary/FiniteElement.o \
	src/ElementLibrary/SpecializedFormulations/UniformStrainHexahedron.o \
	src/ElementLibrary/SpecializedFormulations/UniformStrainQuadrilateral.o \
	src/ElementLibrary/SpecializedFormulations/SimpleTetrahedron.o \
	src/ElementLibrary/SpecializedFormulations/LinearTriangle.o \
	src/ElementLibrary/SpecializedFormulations/TriangleShell.o \
	src/ElementLibrary/SpecializedFormulations/Line.o \
	src/ElementLibrary/SpecializedFormulations/QuadrilateralShell.o \
	src/ElementLibrary/FiniteElementUtilities.o

IC_OBJS ?=\
	src/InitialConditions/InitialConditions.o

IO_OBJS ?=\
	src/IO/LogStream.o \
	src/IO/BinStream.o \
	src/IO/FileManagerDataT.o \
	src/IO/AbaqusFileManagerDataT.o \
	src/IO/EllipsoidFileManagerDataT.o \
	src/IO/FileManagerT.o \
	src/IO/PMPIOBase.o \
	src/IO/RestartFile.o \
	src/IO/silo/SiloFile.o \
	src/IO/TimeHistoryFileManager.o \
	src/IO/ticpp/HierarchicalDataNode.o \
	src/IO/ticpp/tinyxml.o \
	src/IO/ticpp/tinyxmlparser.o \
	src/IO/ticpp/tinyxmlerror.o 

LINALG_OBJS ?=\
	src/LinearAlgebra/TrilinosVector.o \
	src/LinearAlgebra/TrilinosDoFHandler.o

# Leave blank. Used to add local files to build during development.
LOCAL_OBJS ?= src/PhysicsSolvers/SW4Solver.o \
              src/PhysicsSolvers/SW4/Source.o \
              src/PhysicsSolvers/SW4/Filter.o \
              src/PhysicsSolvers/SW4/GridPointSource.o \
              src/PhysicsSolvers/SW4/Qspline.o \
              src/PhysicsSolvers/SW4/SuperGrid.o \
              src/PhysicsSolvers/SW4/time_functions.o \
              src/PhysicsSolvers/SW4/MaterialBlock.o \
              src/PhysicsSolvers/SW4/TimeSeries.o \
              src/PhysicsSolvers/SW4/MaterialPfile.o

MATERIAL_OBJS ?=\
	src/Constitutive/ConstitutivePropertiesTable.o \
	src/Constitutive/Material/MaterialFactory.o \
	src/Constitutive/Material/MaterialBase.o \
	src/Constitutive/Material/LinearElasticDEM.o \
	src/Constitutive/Material/LinearElasticIntermediate.o \
	src/Constitutive/Material/LinearElastic.o \
	src/Constitutive/CohesiveZone/CohesiveZoneBase.o \
	src/Constitutive/CohesiveZone/InitiallyRigidCohesiveZone.o \
	src/Constitutive/CohesiveZone/CohesiveZoneFactory.o \
	src/Constitutive/Interface/InterfaceFactory.o \
	src/Constitutive/Interface/InterfaceBase.o \
	src/Constitutive/Interface/HertzianIntermediate.o \
	src/Constitutive/Interface/Hertzian.o \
	src/Constitutive/Interface/Linearized.o \
	src/Constitutive/Interface/PenaltyCoulombIntermediate.o \
	src/Constitutive/Interface/PenaltyCoulomb.o \

GMM_OBJS =

MISC_OBJS ?=\
	src/Common/Common.o \
	src/Common/Version.o \
	src/Main/main.o \
	src/SurfaceGeneration/FractunatorFactory.o \
	src/SurfaceGeneration/FractunatorBase.o \
	src/SurfaceGeneration/JointSetT.o \
	src/SurfaceGeneration/JointPopulator.o \
	src/SurfaceGeneration/Geostatistics.o \
	src/SurfaceGeneration/StatisticalDistributionBaseT.o \
	src/Interpolation/SurfaceKernel.o \
	src/Interpolation/VolumeKernel.o \
	src/SurfaceGeneration/FractalBaseT.o \
	src/SurfaceGeneration/FractalSurface.o \
	src/SurfaceGeneration/FractalVolume.o

MPICOMM_OBJS ?=\
	src/MPI_Communications/NeighborCommunication.o \
	src/MPI_Communications/PartitionBase.o \
	src/MPI_Communications/SpatialPartition.o \
	src/MPI_Communications/Bifroest/BifroestBase.o

OBJMANAGER_OBJS ?=\
	src/ObjectManagers/PhysicalDomainT.o \
	src/ObjectManagers/ElementManagerT.o \
	src/ObjectManagers/ProblemManagerT.o \
	src/ObjectManagers/EdgeManagerT.o \
	src/ObjectManagers/FaceManagerT.o \
	src/ObjectManagers/ExternalFaceManagerT.o \
	src/ObjectManagers/DiscreteElementManagerBaseT.o \
	src/ObjectManagers/DiscreteElementManagerT.o \
	src/ObjectManagers/EllipsoidalDiscreteElementManagerT.o \
	src/ObjectManagers/ContactManagerBaseT.o \
	src/ObjectManagers/EllipsoidalContactManagerT.o \
	src/ObjectManagers/ContactManagerT.o \
	src/ObjectManagers/TableManager.o \
	src/ObjectManagers/CartesianGridManagerT.o \
	src/ObjectManagers/EnergyT.o \
	src/ObjectManagers/UnitManager.o

SOLVERS_OBJS ?=\
	src/PhysicsSolvers/SolverApplicationSet.o \
	src/PhysicsSolvers/SolverBase.o \
	src/PhysicsSolvers/SolverFactory.o \
	src/PhysicsSolvers/ElementTester.o \
	src/PhysicsSolvers/UpdateFieldWithFunction.o  \
	src/PhysicsSolvers/UpdateParentIndicies.o  \
	src/PhysicsSolvers/ParallelPlateFlowSolverBase.o \
	src/PhysicsSolvers/ParallelPlateFlowSolver.o \
	src/PhysicsSolvers/HydroStaticParallelPlateFlowSolver.o \
	src/PhysicsSolvers/Flow/ParallelPlateFlowSolverFV.o  \
	src/PhysicsSolvers/Lagrange/LagrangeSolverBase.o \
	src/PhysicsSolvers/Lagrange/LagrangeLargeStrain.o \
	src/PhysicsSolvers/Lagrange/LagrangeSmallStrain.o \
	src/PhysicsSolvers/Lagrange/LagrangeSmallStrainLinearElastic.o \
	src/PhysicsSolvers/Lagrange/LagrangeHelperFunctions.o \
	src/PhysicsSolvers/LagrangeDynamicsParallelPlateFlowExplicit.o \
	src/PhysicsSolvers/LagrangeExplicitDynamicsSolver.o \
	src/PhysicsSolvers/ImplicitLaplace.o \
	src/PhysicsSolvers/ImplicitMechanics.o \
	src/PhysicsSolvers/SubstepSolver.o \
	src/PhysicsSolvers/WriteFieldToFile.o \
	src/PhysicsSolvers/Lagrange/MechanicsKernel.o
	
TOTALVIEW_OBJS =

UTILITIES_OBJS ?=\
	src/Utilities/Kinematics.o \
	src/Utilities/Functions.o \
	src/Utilities/StringUtilities.o \
	src/Utilities/RCVSparse.o  \
	src/Utilities/TrilinosUtilities.o \
	src/Utilities/GeometryUtilities.o \
	src/MeshUtilities/MeshGenerator.o \
	src/MeshUtilities/SimpleGeometricObjects.o \
	src/MeshUtilities/MeshUtilities.o


############ Uniaxial extension using transversely isotropic cube ######################
## This script solves a uniaxial extension on a cubic which is transversely isotropic
## with homogeneous fibre directions. 
## The geometry of the cube is interpolated using a tricubic hermite element. 


### Step 0: Housekeeping ##############################################################
import os, sys
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')));
from math import pi
from opencmiss import CMISS

# Problem parameters
length = 1.0
width = 1.0
height = 1.0

numGlobalXElements = 1
numGlobalYElements = 1
numGlobalZElements = 1
numOfXi = 3

# User numbers
(coordinateSystemUserNumber,
	regionUserNumber,
	linearBasisUserNumber,
	cubicBasisUserNumber,
	generatedMeshUserNumber,
	meshUserNumber,
	decompositionUserNumber,
	geometricFieldUserNumber,
	fibreFieldUserNumber,	
	materialFieldUserNumber,
	dependentFieldUserNumber,
	deformedFieldUserNumber,
	equationsSetFieldUserNumber,
	equationsSetUserNumber,
	equationsUserNumber,
	problemUserNumber) = range(1,17)

### Step 1: Set up parallel computing ##################################################
numberOfNodes = CMISS.ComputationalNumberOfNodesGet()
rankNumber = CMISS.ComputationalNodeNumberGet()


### Step 1: Coordinate system ##########################################################
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

### Step 2: Region #####################################################################
region = CMISS.Region()
region.CreateStart(regionUserNumber, CMISS.WorldRegion)
region.CoordinateSystemSet(coordinateSystem)
region.LabelSet("Region")
region.CreateFinish()

### Step 3: Basis functions ############################################################
linearBasis = CMISS.Basis()
linearBasis.CreateStart(linearBasisUserNumber)
linearBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numOfXi)
linearBasis.QuadratureNumberOfGaussXiSet([CMISS.BasisQuadratureSchemes.HIGH]*numOfXi)
linearBasis.QuadratureLocalFaceGaussEvaluateSet(True)
linearBasis.CreateFinish()

cubicBasis = CMISS.Basis()
cubicBasis.CreateStart(cubicBasisUserNumber)
cubicBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE]*numOfXi)
cubicBasis.QuadratureNumberOfGaussXiSet([CMISS.BasisQuadratureSchemes.HIGH]*numOfXi)
cubicBasis.QuadratureLocalFaceGaussEvaluateSet(True)
cubicBasis.CreateFinish()

cubicBasisNumber = 1
linearBasisNumber = 2
bases = [cubicBasis, linearBasis]

### Step 5: Set up mesh ###############################################################
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.TypeSet(CMISS.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet(bases)
generatedMesh.ExtentSet([width, length, height])
generatedMesh.NumberOfElementsSet([numGlobalXElements, numGlobalYElements, numGlobalZElements])
mesh = CMISS.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

### Step 6: Decomposition #############################################################
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfNodes)
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

### Step 7: Create geometric field ####################################################
geometricField = CMISS.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Geometry")
geometricField.NumberOfVariablesSet(1)
geometricField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)
for i in range (1,4):
	geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicBasisNumber)
geometricField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
geometricField.CreateFinish()
# Update geometric parameters
generatedMesh.GeometricParametersCalculate(geometricField)

### Step 8: Create fibre field #######################################################
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre")
fibreField.NumberOfVariablesSet(1)
fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)
fibreField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
#for component in range (1,4):
	#fibreField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, component, CMISS.FieldInterpolationTypes.CONSTANT)
fibreField.CreateFinish()

# For future reference, how to get the nodes in a particular component of mesh. 
#nodes = CMISS.Nodes()
#region.NodesGet(nodes)
#numberOfNodes = nodes.NumberOfNodesGet()
#cubicNodes = [n for n in range(1,numberOfNodes + 1) if mesh.NodeExists(cubicBasisNumber,n)]

# Initialise the fibre rotation angles in radians
totalNodes = [1,2,3,4,5,6,7,8]
fibreFieldAngle = 30 * pi /180  # 30 degrees in radians
fibreAngles = [fibreFieldAngle, 0, 0]
for node in totalNodes:
	for component, fibreAngle in enumerate(fibreAngles, 1):
		fibreField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, node, component, fibreAngle)		
		#fibreField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, fibreAngle)

### Step 9: Material field ##########################################################
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber, region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
materialField.MeshDecompositionSet(decomposition)
materialField.GeometricFieldSet(geometricField)
materialField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material")
materialField.NumberOfVariablesSet(1)
materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 5)
for i in range (1,5):
	materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicBasisNumber)
materialField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
materialField.CreateFinish()
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 2.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2, 5.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 3, 10.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, 0.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 5, 5.0)

### Step 10: Dependent field ########################################################
dependentField = CMISS.Field()
dependentField.CreateStart(dependentFieldUserNumber, region)
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Dependent")
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT)
dependentField.NumberOfVariablesSet(2)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN, 4)
for i in range (1,4):
	dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicBasisNumber)
	dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i, cubicBasisNumber)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 4, linearBasisNumber)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 4, linearBasisNumber)
dependentField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry. 
for i in range (1,4):
	CMISS.Field.ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, i, 
dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, i)

# Set hydrostatic pressure. 
dependentField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -8.0)

### Step 11: Create Equation set #####################################################
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_EXPONENTIAL, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(CMISS.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

### Step 12: Problem ################################################################
problem =CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY, CMISS.ProblemTypes.FINITE_ELASTICITY, CMISS.ProblemSubTypes.NONE)
problem.CreateFinish()

problem.ControlLoopCreateStart()
controlLoop = CMISS.ControlLoop()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
problem.ControlLoopCreateFinish()

### Step 13: Solvers ################################################################
nonLinearSolver = CMISS.Solver()
linearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.OutputTypeSet(CMISS.SolverOutputTypes.PROGRESS)
nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.EQUATIONS)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.LinearTypeSet(CMISS.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

### Step 14: Solver equations #######################################################
solver = CMISS.Solver()
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverEquationsGet([CMISS.ControlLoopIdentifiers.NODE], 1, solverEquations)
solverEquations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

### Step 15: Boundary conditions #####################################################
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

# Set face 1,3,5,7 fixed in all directions.
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,3,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,7,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,5,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,7,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)

# Set displacement of face 2,4,6,8 in x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,2,1,CMISS.BoundaryConditionsTypes.FIXED,-5.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,4,1,CMISS.BoundaryConditionsTypes.FIXED,-5.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,6,1,CMISS.BoundaryConditionsTypes.FIXED,-5.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,8,1,CMISS.BoundaryConditionsTypes.FIXED,-5.0)

# Set single arc-length derivatives
for node in range (1,9):
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,node,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,node,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,node,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,node,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,node,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,node,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)

# Set all cross derivatives to zero in all three directions. 
for component in range (1,4):
	for node in range (1,9):
		boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,node,component,CMISS.BoundaryConditionsTypes.FIXED,0.0)
		boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,node,component,CMISS.BoundaryConditionsTypes.FIXED,0.0)
		boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,node,component,CMISS.BoundaryConditionsTypes.FIXED,0.0)
		boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,node,component,CMISS.BoundaryConditionsTypes.FIXED,0.0)

solverEquations.BoundaryConditionsCreateFinish()


### Step 14: Solve the problem ########################################################
problem.Solve()

### Step 15: Housekeeping #############################################################
# Copy deformed geometry into a new geometric field for visualisation. 
deformedField = CMISS.Field()
deformedField.CreateStart(deformedFieldUserNumber, region)
deformedField.MeshDecompositionSet(decomposition)
deformedField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
deformedField.VariableLabelSet(CMISS.FieldVariableTypes.U, "DeformedGeometry")
for component in [1,2,3,]:
	deformedField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, component, cubicBasisNumber)
deformedField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
deformedField.CreateFinish()
for component in [1,2,3]:
	dependentField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, deformedField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)


# Export results to .exnode and .exelem files. 
exportFields = CMISS.Fields()
exportFields.CreateRegion(region)
exportFields.NodesExport("../Results/UniAxialTransverselyIsotropic","FORTRAN")
exportFields.ElementsExport("../Results/UniAxialTransverselyIsotropic","FORTRAN")
exportFields.Finalise()



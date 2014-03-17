############ Uniaxial extension using transversely isotropic cube ######################
## This script solves a uniaxial extension on a cubic which is transversely isotropic
## with homogeneous fibre directions. 
## The geometry of the cube is interpolated using a tricubic hermite element. 


### Step 0: Housekeeping ##############################################################
import os, sys
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')));
from opencmiss import CMISS
CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])


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
	quadraticBasisUserNumber,
	cubicBasisUserNumber,
	generatedMeshUserNumber,
	meshUserNumber,
	decompositionUserNumber,
	geometricFieldUserNumber,
	fibreFieldUserNumber,	
	materialFieldUserNumber,
	dependentFieldUserNumber,
	equationsSetFieldUserNumber,
	equationsSetUserNumber,
	equationsUserNumber,
	problemUserNumber) = range(1,16)

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
quadraticBasis = CMISS.Basis()
quadraticBasis.CreateStart(quadraticBasisUserNumber)
quadraticBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numOfXi)
quadraticBasis.QuadratureNumberOfGaussXiSet([CMISS.BasisQuadratureSchemes.HIGH]*numOfXi)
quadraticBasis.QuadratureLocalFaceGaussEvaluateSet(True)
quadraticBasis.CreateFinish()

cubicBasis = CMISS.Basis()
cubicBasis.CreateStart(cubicBasisUserNumber)
cubicBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE]*numOfXi)
cubicBasis.QuadratureNumberOfGaussXiSet([CMISS.BasisQuadratureSchemes.HIGH]*numOfXi)
cubicBasis.QuadratureLocalFaceGaussEvaluateSet(True)
cubicBasis.CreateFinish()
quadraticBasisNumber = 1
cubicBasisNumber = 2
bases = [quadraticBasis, cubicBasis]

### Step 5: Set up mesh ###############################################################
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.TypeSet(CMISS.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet(bases) 
generatedMesh.ExtentSet([width, length, height])
generatedMesh.NumberOfElementsSet([numGlobalXElements,numGlobalYElements,numGlobalZElements]) 
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
geometricField.CreateFinish()
# Update geometric parameters
generatedMesh.GeometricParametersCalculate(geometricField)


### Step 8: Create fibre field #######################################################
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreFieldNumOfComponents = 1
fibreField.NumberOfVariablesSet(1)
fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, fibreFieldNumOfComponents)
for i in range (1,fibreFieldNumOfComponents +1):
	fibreField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i,cubicBasisNumber)
fibreField.CreateFinish()
# Initialise the fibre rotation angles
fibreFieldAngle = 0.523598776  # 30 degrees in radians
for node in range (1,9):
	for component in range(1,fibreFieldNumOfComponents + 1):
		fibreField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, node, component, fibreFieldAngle)

### Step 9: Material field ##########################################################
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber, region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
materialField.GeometricFieldSet(geometricField)
materialField.NumberOfVariablesSet(1)
materialField.NumberOfComponentsSet(5)
for i in range (1,6):
	materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicBasisNumber)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, 1, 1.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, 2, 5.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, 3, 10.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, 4, 0.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, 5, 5.0)

### Step 10: Dependent field ########################################################
dependentField = CMISS.Field()
dependentField.CreateStart(dependentFieldUserNumber, region)
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField)
dependentField.DependentTypeSet(CMISS.DependentFieldTypes.DEPENDENT)
dependentField.NumberOfVariablesSet(2)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN, 4)
for i in range (1,4):
	dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicBasisNumber)
	dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i, cubicBasisNumber)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 4, quadraticBasisNumber)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 4, quadraticBasisNumber)
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry. 
for i in range (1,4):
	ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, i, 
dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, i)
# Set hydrostatic pressure. 
dependentField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldSetTypes.VALUES, 4, -14.0)

### Step 11: Create Equation set #####################################################
equationsSetField = CMISS.Field()
equationsSetField.CreateStart(equationsSetFieldUserNumber, region)
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, FibreField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_GUCCIONE, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(CMISS.EquationsOutputTypes.NO_OUTPUT)
equationsSet.EquationsCreateFinish()

### Step 12: Problem ################################################################
problem =CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY, CMISS.ProblemTypes.FINITE_ELASTICITY, CMISS.ProblemSubTypes.NO_SUBTYPE)
problem.CreateFinish()

problem.ControlLoopCreateStart()
controlLoop = CMISS.ControlLoop()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
#controlLoop.MaximumIterationsSet(100)
problem.ControlLoopCreateFinish()

### Step 13: Solvers ################################################################
nonLinearSolver = CMISS.Solver()
linearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet(CMISS.ControlLoopIdentifiers.NODE, 1, nonLinearSolver)
nonLinearSolver.OutputTypeSet(CMISS.SolverOutput.PROGRESS)
nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianEquationsTypes.EQUATIONS)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.TypeSet(CMISS.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

### Step 14: Solver equations #######################################################
solver = CMISS.Solver()
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverEquationsGet(solverEquations)
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
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,2,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,4,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,6,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,8,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)

# Set single arc-length derivatives
# Node 1
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,1,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 2
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,2,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,2,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 3
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,3,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,3,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 4
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,4,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,4,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 5
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,5,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,5,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 6
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,6,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,6,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 7
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,7,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,7,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 8
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,8,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3,8,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)

# Set all cross derivatives to zero in x direction. 
# Node 1
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 2
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 3
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 4
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 5
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 6
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 7
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 8
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)

# Set all cross derivatives to zero in y direction.
# Node 1
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 2
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 3
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 4
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 5
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 6
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 7
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 8
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)

# Set all cross derivatives to zero in z direction.
# Node 1
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 2
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 3
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 4
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 5
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 6
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 7
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 8
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)

solverEquations.BoundaryConditionsCreateFinish()


### Step 14: Solve the problem ########################################################
problem.Solve()

### Step 15: Housekeeping #############################################################
exportFields = CMISS.Fields()
exportFields.CreateRegion(region)
exportFields.NodesExport("UniAxialTransverselyIsotropic","FORTRAN")
exportFields.ElementsExport("UniAxialTransverselyIsotropic","FORTRAN")
exportFields.Finalise()

######## Uniaxial extension using tricubic hermite element ###############
## This script file solves a uniaxial extension in the x axis on a single
## tricubic hermite cuboid element. 

### Step 0: Housekeeping
import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')));

# Initialise OpenCMISS
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
	meshUserNumber,
	decompositionUserNumber,
	geometricFieldUserNumber,
	fibreFieldUserNumber,	
	materialFieldUserNumber,
	dependentFieldUserNumber,
	equationsSetFieldUserNumber,
	equationsSetUserNumber,
	equationsUserNumber,
	problemUserNumber) = range(1,15)


### Step 1: Allow for parallel computation ############################################
numberOfNodes = CMISS.ComputationalNumberOfNodesGet()
rankNumber = CMISS.ComputationalNodeNumberGet()

### Step 2: Set up co-ordinate system #################################################
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

### Step 3: Set up region #############################################################
region = CMISS.Region()
region.CreateStart(regionUserNumber, CMISS.WorldRegion)
region.CoordinateSystemSet(coordinateSystem)
region.LabelSet("Region")
region.CreateFinish()

### Step 4: Set up tricubic hermite basis functions ###################################
linearBasis = CMISS.Basis()  # For interpolation of hydrostatic pressure. 
linearBasis.CreateStart(linearBasisUserNumber)
#linearBasis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
#linearBasis.NumberOfXiSet(numOfXi)
#linearBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE]*3)
linearBasis.QuadratureNumberOfGaussXiSet([3]*numOfXi)
linearBasis.CreateFinish()

cubicBasis = CMISS.Basis() # For geometry. 
cubicBasis.CreateStart(cubicBasisUserNumber)
cubicBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE]*3)
cubicBasis.QuadratureNumberOfGaussXiSet([3]*numOfXi)
cubicBasis.QuadratureLocalFaceGaussEvaluateSet(True)
cubicBasis.CreateFinish()

### Step 5: Set up mesh ###############################################################
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber, region, numOfXi)
mesh.NumberOfComponentsSet(2)
mesh.NumberOfElementsSet(1)

nodes = CMISS.Nodes()
nodes.CreateStart(region, 8)
nodes.CreateFinish()

# Linear lagrange component
linearMeshComponentNumber = 1
linearElem = CMISS.MeshElements()
linearElem.CreateStart(mesh,linearMeshComponentNumber,linearBasis)
linearElem.NodesSet(1,[1,2,3,4,5,6,7,8])
linearElem.CreateFinish()

# Cubic Hermite component
cubicMeshComponentNumber = 2
cubicElem = CMISS.MeshElements()
cubicElem.CreateStart(mesh,cubicMeshComponentNumber,cubicBasis)
cubicElem.NodesSet(1,[1,2,3,4,5,6,7,8])
cubicElem.CreateFinish()
mesh.CreateFinish()

### Step 6: Create decomposition ######################################################
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfNodes)
decomposition.CreateFinish()

### Step 7: Create geometric field ####################################################
geometricField = CMISS.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
geometricField.MeshDecompositionSet(decomposition)
geometricField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Geometry")
for i in range(1,4):
	geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicMeshComponentNumber)
geometricField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
geometricField.CreateFinish()
# Node 1
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,1,1,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,1,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,1,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,1,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,1,1,3,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,1,3,1.0)

# Node 2
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,2,1,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,2,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,2,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,2,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,2,3,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,2,3,1.0)

# Node 3
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,3,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,3,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,3,2,length)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,3,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,3,3,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,3,3,1.0)
# Node 4
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,4,1,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,4,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,4,2,length)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,4,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,4,3,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,4,3,1.0)

# Node 5
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,5,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,5,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,5,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,5,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,5,3,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,5,3,1.0)

# Node 6
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,6,1,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,6,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,6,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,6,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,6,3,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,6,3,1.0)
# Node 7
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,7,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,7,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,7,2,length)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,7,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,7,3,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,7,3,1.0)

# Node 8
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,8,1,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2,8,1,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,8,2,length)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,3,8,2,1.0)

geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,1,8,3,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,5,8,3,1.0)

# Export geometric field for debugging. 
exportGeometricFields = CMISS.Fields()
exportGeometricFields.CreateRegion(region)
exportGeometricFields.NodesExport("GeometricField","FORTRAN")
exportGeometricFields.ElementsExport("GeometricField","FORTRAN")
exportGeometricFields.Finalise()

### Step 8: Create fibre field ########################################################
fibre = CMISS.Field()
fibre.CreateStart(fibreFieldUserNumber, region)
fibre.TypeSet(CMISS.FieldTypes.FIBRE)
fibre.MeshDecompositionSet(decomposition)
fibre.GeometricFieldSet(geometricField)
fibre.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre")
for i in range(1,4):
	fibre.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, linearMeshComponentNumber)
fibre.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
fibre.CreateFinish()

### Step 9: Create material field #####################################################
material = CMISS.Field()
material.CreateStart(materialFieldUserNumber, region)
material.TypeSet(CMISS.FieldTypes.MATERIAL)
material.MeshDecompositionSet(decomposition)
material.GeometricFieldSet(geometricField)
material.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material")
material.NumberOfVariablesSet(1)
material.NumberOfComponentsSet(CMISS.FieldVariableTypes.U,2)
material.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 1, linearMeshComponentNumber)
material.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 2, linearMeshComponentNumber)
material.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
material.CreateFinish() 

material.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,1,2.0)
material.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,2,4.0)


### Step 10: Create dependent geometric field #########################################
dependentField = CMISS.Field()
dependentField.CreateStart(dependentFieldUserNumber, region)
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Dependent")
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT)
dependentField.NumberOfVariablesSet(2)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN,4)
for i in range(1,4):
	dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i,cubicMeshComponentNumber)
	dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i,cubicMeshComponentNumber)

dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 4, linearMeshComponentNumber)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 4, linearMeshComponentNumber)
dependentField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
dependentField.CreateFinish()
# Initialise with undeformed geometry

for i in range(1,4):
	CMISS.Field.ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, i, dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, i)

dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -8.0)

### Step 11: Create Equations Set and equations #######################################
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.MOONEY_RIVLIN, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, material)
equationsSet.MaterialsCreateFinish()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(CMISS.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

### Step 12: Create Problem ###########################################################
# Define problem 
problem = CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY, CMISS.ProblemTypes.FINITE_ELASTICITY, CMISS.ProblemSubTypes.NONE)
problem.CreateFinish()
# Create control loops
problem.ControlLoopCreateStart()
controlLoop = CMISS.ControlLoop()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.MaximumIterationsSet(2)
problem.ControlLoopCreateFinish()
# Create nonlinear numerical solver
linearSolver = CMISS.Solver()
nonLinearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.OutputTypeSet(CMISS.SolverOutputTypes.PROGRESS)
nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonRelativeToleranceSet(1.0E-10)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.LinearTypeSet(CMISS.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()
# Add solver equations sets which encompass the physics
solverEquations = CMISS.SolverEquations()
solver = CMISS.Solver()
problem.SolverEquationsCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(CMISS.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

### Step 13: Prescribe Boundary Conditions ############################################
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
# Set face 1,3,5,7 fixed in all directions.
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)

boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,3,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,7,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)

boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,5,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,7,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)

# Set displacement of face 2,4,6,8 in x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,2,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,4,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,6,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,8,1,CMISS.BoundaryConditionsTypes.FIXED,1.1)

# Set single arc-length derivatives
# Node 1
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,1,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,1,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 2
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,2,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,2,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 3
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,3,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,3,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 4
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,4,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,4,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 5
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,5,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,5,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 6
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,6,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,6,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 7
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,7,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,7,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)
# Node 8
# x direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# y direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,8,2,CMISS.BoundaryConditionsTypes.FIXED,1.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# z direction
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,3,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,5,8,3,CMISS.BoundaryConditionsTypes.FIXED,1.0)

# Set all cross derivatives to zero in x direction. 
# Node 1
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 2
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 3
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 4
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 5
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 6
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 7
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 8
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)

# Set all cross derivatives to zero in y direction.
# Node 1
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 2
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 3
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,3,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 4
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,4,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 5
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 6
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 7
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,7,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 8
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,8,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)

# Set all cross derivatives to zero in z direction.
# Node 1
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 2
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 3
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 4
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 5
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,5,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 6
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,6,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 7
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,7,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
# Node 8
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,4,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,6,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,7,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1,8,8,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)

solverEquations.BoundaryConditionsCreateFinish()

### Step 14: Solve the problem ########################################################
problem.Solve()

### Step 15: Housekeeping #############################################################
exportFields = CMISS.Fields()
exportFields.CreateRegion(region)
exportFields.NodesExport("UniaxialExtensionTricubicHermite","FORTRAN")
exportFields.ElementsExport("UniaxialExtensionTricubicHermite","FORTRAN")
exportFields.Finalise()

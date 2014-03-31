######## Uniaxial extension using tricubic hermite element ###############
## This script file solves a uniaxial extension in the x axis on a single
## tricubic hermite cuboid element. 

### Step 0: Housekeeping
import sys, os
from math import pi
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
	deformedFieldUserNumber, 
	dependentFieldUserNumber,
	equationsSetFieldUserNumber,
	equationsSetUserNumber,
	equationsUserNumber,
	problemUserNumber) = range(1,16)


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
linearBasis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
linearBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3)
linearBasis.QuadratureNumberOfGaussXiSet([CMISS.BasisQuadratureSchemes.MID]*numOfXi)
linearBasis.CreateFinish()

cubicBasis = CMISS.Basis() # For geometry. 
cubicBasis.CreateStart(cubicBasisUserNumber)
cubicBasis.TypeSet(CMISS.BasisTypes.LAGRANGE_HERMITE_TP)
cubicBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE]*3)
cubicBasis.QuadratureNumberOfGaussXiSet([CMISS.BasisQuadratureSchemes.MID]*numOfXi)
cubicBasis.QuadratureLocalFaceGaussEvaluateSet(True)
cubicBasis.CreateFinish()

cubicMeshComponentNumber = 1
linearMeshComponentNumber = 2
### Step 5: Set up mesh ###############################################################
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber, region, numOfXi)
mesh.NumberOfComponentsSet(2)
mesh.NumberOfElementsSet(1)

nodes = CMISS.Nodes()
nodes.CreateStart(region, 8)
nodes.CreateFinish()
# Cubic Hermite component
cubicElem = CMISS.MeshElements()
cubicElem.CreateStart(mesh,cubicMeshComponentNumber,cubicBasis)
cubicElem.NodesSet(1,[1,2,3,4,5,6,7,8])
cubicElem.CreateFinish()

# Linear lagrange component
linearElem = CMISS.MeshElements()
linearElem.CreateStart(mesh,linearMeshComponentNumber,linearBasis)
linearElem.NodesSet(1,[1,2,3,4,5,6,7,8])
linearElem.CreateFinish()
mesh.CreateFinish()
### Step 6: Create decomposition ######################################################
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfNodes)
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
    geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicMeshComponentNumber)
geometricField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
geometricField.CreateFinish()

# Initialise geometric parameters
xNodes = [0.0, 1.0, 0.0, 1.0, 0.0,  1.0, 0.0, 1.0]
yNodes = [0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0]
zNodes = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]

for node, value in enumerate(xNodes, 1):
    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, value)
for node, value in enumerate(yNodes, 1):
    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2, value)
for node, value in enumerate(zNodes, 1):
    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3, value)

# Set all single derivatives 
Deriv = [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
for node in range (1,9):
    for component in range (1,4):
	    for deriv in Deriv:
		    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, deriv, node, component, 0.0)

for node in range (1,9):
    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1, node, 1, 1.0)
    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, node, 3, 1.0)
    geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3, node, 2, -1.0)

# Export geometric field for debugging. 
exportGeometricFields = CMISS.Fields()
exportGeometricFields.CreateRegion(region)
exportGeometricFields.NodesExport("GeometricField","FORTRAN")
exportGeometricFields.ElementsExport("GeometricField","FORTRAN")
exportGeometricFields.Finalise()

### Step 8: Create fibre field ########################################################
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre")
fibreField.NumberOfVariablesSet(1)
fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)
for i in range (1,4):
	fibreField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, linearMeshComponentNumber)
fibreField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
fibreField.CreateFinish()

# Initialise the fibre rotation angles in radians
fibreAngle = [0,0,0]
for component, fibre in enumerate(fibreAngle,1):
    fibreField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, fibre)

### Step 9: Create material field #####################################################
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber, region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
materialField.MeshDecompositionSet(decomposition)
materialField.GeometricFieldSet(geometricField)
materialField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material")
materialField.NumberOfVariablesSet(1)
materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 5)
for i in range (1,6):
    materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, linearMeshComponentNumber)
materialField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
materialField.CreateFinish()

materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 2.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2, 5.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 3, 10.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, 0.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 5, 5.0)

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

for component in range(1,4):
    CMISS.Field.ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)

dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -10.0)

### Step 11: Create Equations Set and equations #######################################
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
problem.ControlLoopCreateFinish()
# Create nonlinear numerical solver
linearSolver = CMISS.Solver()
nonLinearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.OutputTypeSet(CMISS.SolverOutputTypes.PROGRESS)
nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonRelativeToleranceSet(1.0E-5)
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
leftFaceNodes = [1,3,5,7]
rightFaceNodes = [2,4,6,8]
bottomFaceNodes = [1,2,5,6]
backFaceNodes = [1,2,3,4]

# Set right face with force application in x direction
for node in rightFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, CMISS.BoundaryConditionsTypes.FIXED, -2.0)

# Set left face fixed in x direction
for node in leftFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set back face fixed in y direction. 
for node in backFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set bottom face fixed in z direction. 
for node in bottomFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U,1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set single derivatives in back face
for node in backFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 1.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set single derivatives in left face
for node in leftFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3, node, 2, CMISS.BoundaryConditionsTypes.FIXED, -1.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 1.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set single derivatives in bottom face
for node in bottomFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3, node, 2, CMISS.BoundaryConditionsTypes.FIXED, -1.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    #boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set all cross derivatives to zero. 
Deriv = [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]  
for node in range (1,9):
    for component in [1,2,3]:
    	for deriv in Deriv:
		boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, deriv, node, component, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

solverEquations.BoundaryConditionsCreateFinish()

### Step 14: Solve the problem ########################################################
problem.Solve()

### Step 15: Housekeeping #############################################################
deformedField = CMISS.Field()
deformedField.CreateStart(deformedFieldUserNumber, region)
deformedField.MeshDecompositionSet(decomposition)
deformedField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
deformedField.VariableLabelSet(CMISS.FieldVariableTypes.U, "DeformedGeometry")
for component in [1,2,3]:
    deformedField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, component, 1)
deformedField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
deformedField.CreateFinish()

for component in [1,2,3]:
    dependentField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, deformedField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)

exportFields = CMISS.Fields()
exportFields.CreateRegion(region)
exportFields.NodesExport("../Results/UniaxialExtensionTricubicHermite_FreeBackFace", "FORTRAN")
exportFields.ElementsExport("../Results/UniaxialExtensionTricubicHermite_FreeBackFace", "FORTRAN")
exportFields.Finalise()

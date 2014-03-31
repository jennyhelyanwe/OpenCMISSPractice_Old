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
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber, region, 3)
mesh.NumberOfComponentsSet(2)
mesh.NumberOfElementsSet(1)

nodes = CMISS.Nodes()
nodes.CreateStart(region, 8)
nodes.CreateFinish()

elemC = CMISS.MeshElements()
elemC.CreateStart(mesh, 1, cubicBasis)
elemC.NodesSet(1,[1,2,3,4,5,6,7,8])
elemC.CreateFinish()

elemL = CMISS.MeshElements()
elemL.CreateStart(mesh, 2, linearBasis)
elemL.NodesSet(1,[1,2,3,4,5,6,7,8])
elemL.CreateFinish()

mesh.CreateFinish()

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
# Set all single derivatives in all components to be 1.0. 
for deriv in [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2, CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
    for node in range (1,9):
    	for component in [1,2,3]:
    		geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, deriv, node,component, 0.0)

# Export geometric field for debugging. 
exportGeometricFields = CMISS.Fields()
exportGeometricFields.CreateRegion(region)
exportGeometricFields.NodesExport("GeometricField","FORTRAN")
exportGeometricFields.ElementsExport("GeometricField","FORTRAN")
exportGeometricFields.Finalise()

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

fibreField.CreateFinish()

# Initialise the fibre rotation angles in radians
totalNodes = [1,2,3,4,5,6,7,8]
fibreAngles = [0, 0, 0]
for node in totalNodes:
	for component, fibreAngle in enumerate(fibreAngles, 1):
		fibreField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, component, fibreAngle)		
		

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

leftFaceNodes = [1,3,5,7]
rightFaceNodes = [2,4,6,8]
bottomFaceNodes = [1,2,5,6]
backFaceNodes = [1,2,3,4]

# Set right face with force application in x direction
for node in rightFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,node, 1,CMISS.BoundaryConditionsTypes.FIXED, -2.0)

# Set left face fixed in x direction
for node in leftFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U,1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set bottom face fixed in z direction. 
for node in bottomFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U,1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set back face fixed in y direction. 
for node in backFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U,1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
"""
#Set single derivatives to 
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
		for deriv in [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,  CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]: 
			boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, deriv, node, component, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
"""
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



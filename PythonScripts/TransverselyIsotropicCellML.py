###################### Uniaxial Load on transversely isotropic cube ##############
## This script solves a uniaxial pressure load along the fibre direction of a 
## unit cube. It uses the transversely isotropic constitutive law. 

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
    independentFieldUserNumber, 
	dependentFieldUserNumber,
	deformedFieldUserNumber,
	equationsSetFieldUserNumber,
	equationsSetUserNumber,
	equationsUserNumber,
	problemUserNumber,
    cellMLUserNumber,
    CellMLModelsFieldUserNumber,
    CellMLParametersFieldUserNumber,
    CellMLIntermediateFieldUserNumber) = range(1,22)

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

### Step 5: Set up mesh ###############################################################
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber, region, 3)
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(1)

nodes = CMISS.Nodes()
nodes.CreateStart(region, 8)
nodes.CreateFinish()

elemL = CMISS.MeshElements()
elemL.CreateStart(mesh,1, linearBasis)
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
	geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
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

### Step 8: Fibre field ###############################################################
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre")
fibreField.NumberOfVariablesSet(1)
fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)
for i in range (1,4):
	fibreField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i,1)

fibreField.CreateFinish()

# Initialise the fibre rotation angles in radians
fibreAngle = [0,0,0]
for component, fibre in enumerate(fibreAngle,1):
    fibreField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, fibre)

### Step 9: Material Field ###########################################################

#MaterialVariableTypes = [CMISS.FieldVariableTypes.U, CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldVariableTypes.U1]
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber, region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
materialField.MeshDecompositionSet(decomposition)
materialField.GeometricFieldSet(geometricField)
materialField.NumberOfVariablesSet(1)
#materialField.VariableTypesSet(MaterialVariableTypes)
materialField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material")
materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
for i in range (1,5):
	materialField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, i, CMISS.FieldInterpolationTypes.CONSTANT)
materialField.CreateFinish()
#parameters = [2.0, 3.0]
parameters = [2.0, 3.0, 1.0, 3.0]
for component, param in enumerate(parameters, 1):
    materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, param)

### Step 10: Create equation set ####################################################
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
#equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_GUCCIONE, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.CONSTITUTIVE_LAW_IN_CELLML_EVALUATE, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
equationsSet.MaterialsCreateFinish()
### Step 11: Create dependent field ##############################################
DependentVariableTypes = [CMISS.FieldVariableTypes.U, CMISS.FieldVariableTypes.DELUDELN, CMISS.FieldVariableTypes.U1, CMISS.FieldVariableTypes.U2 ]

dependentField = CMISS.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Dependent")
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT)
dependentField.NumberOfVariablesSet(4)
dependentField.VariableTypesSet(DependentVariableTypes)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U1, 6) # storing strain
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U2, 6) # storing stress
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN, 4)

for component in range (1,7):
    dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1, component, CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2, component, CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)

for i in range (1,5):
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i,1) 
for i in range (1,7):
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U2, i, 1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U1, i, 1)
       
dependentField.CreateFinish()


equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(CMISS.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()
# Initialise from undeformed geometry
for component in [1,2,3]:
    geometricField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,component, dependentField,CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)

dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -8.0)

# Initialise strain and stress fields
for i in range(1,7):
    dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U1, CMISS.FieldParameterSetTypes.VALUES, i, 0.0)
    dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U2, CMISS.FieldParameterSetTypes.VALUES, i, 0.0)

### Step 12: Create CellML environment ###########################################
TransIsoModelIndex = 1
cellML = CMISS.CellML()
cellML.CreateStart(cellMLUserNumber, region)
cellML.ModelImport("transversely_isotropic.cellml")
parameter = ["C1", "C2","C3", "C4"]
strain = ["E11", "E12", "E13", "E22", "E23", "E33"]
stress2PK = ["Tdev11", "Tdev12", "Tdev13", "Tdev22", "Tdev23", "Tdev33"]

# Set strains as known in CellML. These will be fed into the model from iron. 
for i in range(0,6):
    cellML.VariableSetAsKnown(TransIsoModelIndex, "equations/"+strain[i])
for i in range(0,4):
    cellML.VariableSetAsKnown(TransIsoModelIndex, "equations/"+parameter[i])

# Set stresses as unknown in CellML. These will be calculated using the transversely isotropic constitutive model 
for i in range (0,6):
    cellML.VariableSetAsWanted(TransIsoModelIndex, "equations/"+stress2PK[i])
cellML.CreateFinish()

### Step 13: Map the variables to CellML model ###################################
cellML.FieldMapsCreateStart()
# Map the strain from dependentField U1 variable to CellML. 
for i in range(0,6):
    cellML.CreateFieldToCellMLMap(dependentField, CMISS.FieldVariableTypes.U1, i+1, CMISS.FieldParameterSetTypes.VALUES, TransIsoModelIndex, "equations/"+strain[i], CMISS.FieldParameterSetTypes.VALUES)
for i in range(0,4):
    print "equations/"+parameter[i]+"\n"
    cellML.CreateFieldToCellMLMap(materialField, CMISS.FieldVariableTypes.U,i+1, CMISS.FieldParameterSetTypes.VALUES, TransIsoModelIndex, "equations/"+parameter[i], CMISS.FieldParameterSetTypes.VALUES)
# Map the stress from CellML to dependentFieldU2 variable
for i in range(0,6):
    cellML.CreateCellMLToFieldMap(TransIsoModelIndex, "equations/"+stress2PK[i], CMISS.FieldParameterSetTypes.VALUES, dependentField, CMISS.FieldVariableTypes.U2, i+1, CMISS.FieldParameterSetTypes.VALUES)
cellML.FieldMapsCreateFinish()

# Create models field for CellML
CellMLModelsField = CMISS.Field()
cellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber, CellMLModelsField)
cellML.ModelsFieldCreateFinish()

# No need to create a state field since we aren't integrating. 

# Create parameters field for CellML, this is used as the strain field. 
CellMLParametersField = CMISS.Field()
cellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber, CellMLParametersField)
cellML.ParametersFieldCreateFinish()

# Create intermediate field for CellML, this is used as the stress field. 
CellMLIntermediateField = CMISS.Field()
cellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber, CellMLIntermediateField)
cellML.IntermediateFieldCreateFinish()

### Step 14: Problem ################################################################
problem =CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY, CMISS.ProblemTypes.FINITE_ELASTICITY, CMISS.ProblemSubTypes.FINITE_ELASTICITY_CELLML)
problem.CreateFinish()

problem.ControlLoopCreateStart()
controlLoop = CMISS.ControlLoop()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
problem.ControlLoopCreateFinish()

### Step 15: Solvers ################################################################
nonLinearSolver = CMISS.Solver()
linearSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.OutputTypeSet(CMISS.SolverOutputTypes.PROGRESS)
nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.LinearTypeSet(CMISS.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

### Step 16: Solver equations #######################################################
solver = CMISS.Solver()
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverEquationsGet([CMISS.ControlLoopIdentifiers.NODE], 1, solverEquations)
solverEquations.SparsityTypeSet(CMISS.EquationsSparsityTypes.SPARSE)
solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

cellMLSolver = CMISS.Solver()
cellMLEquations = CMISS.CellMLEquations()
problem.CellMLEquationsCreateStart()
nonLinearSolver.NewtonCellMLSolverGet(cellMLSolver)
cellMLSolver.CellMLEquationsGet(cellMLEquations)
cellMLEquations.CellMLAdd(cellML)
problem.CellMLEquationsCreateFinish()

### Step 17: Boundary Conditions ###################################################
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

leftFaceNodes = [1,3,5,7]
rightFaceNodes = [2,4,6,8]
bottomFaceNodes = [1,2,5,6]
backFaceNodes = [1,2,3,4]
# Set left face fixed in x direction
for node in leftFaceNodes:
	boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set right face with force application in x direction

rightFaceNormalXi = 1
BCPressure = CMISS.BoundaryConditionsTypes.PRESSURE_INCREMENTED
for node in rightFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,node,rightFaceNormalXi,BCPressure,-2.097)

# Set right face with x extension
for node in rightFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U, 1, 1,node, 1,CMISS.BoundaryConditionsTypes.FIXED, 1.1)

# Set bottom face fixed in z direction. 
for node in bottomFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U,1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3, CMISS.BoundaryConditionsTypes.FIXED, 0.0)

# Set back face fixed in y direction. 
for node in backFaceNodes:
    boundaryConditions.SetNode(dependentField, CMISS.FieldVariableTypes.U,1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2, CMISS.BoundaryConditionsTypes.FIXED, 0.0)
solverEquations.BoundaryConditionsCreateFinish()

### Step 18: ####################################################################
problem.Solve()

### Step 19: Housekeeping #######################################################
deformedField = CMISS.Field()
deformedField.CreateStart(deformedFieldUserNumber, region)
deformedField.MeshDecompositionSet(decomposition)
deformedField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
deformedField.VariableLabelSet(CMISS.FieldVariableTypes.U, "DeformedGeometry")
for component in [1,2,3]:
	deformedField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, component, 1)
deformedField.CreateFinish()

for component in [1,2,3]:
	dependentField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, deformedField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)

exportFields = CMISS.Fields()
exportFields.CreateRegion(region)
exportFields.NodesExport("../Results/TransverselyIsotropicCellML","FORTRAN")
exportFields.ElementsExport("../Results/TransverselyIsotropicCellML","FORTRAN")
exportFields.Finalise()



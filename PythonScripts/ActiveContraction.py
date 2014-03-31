######################### Active Contraction of Tricubic hermite unit cube ##############
## This script solves an active contraction along the fibre direction of a 
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
	quadraticBasisUserNumber,
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
	problemUserNumber) = range(1,18)

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
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber, region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
materialField.MeshDecompositionSet(decomposition)
materialField.GeometricFieldSet(geometricField)
materialField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material")
materialField.NumberOfVariablesSet(1)
materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 7)
for i in range (1,8):
	materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, cubicBasisNumber)
materialField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
materialField.CreateFinish()

# Set costa material parameters. 
costaParameters = [0.2, 30.0, 12.0, 14.0, 14.0, 10.0, 18.0]
for component, parameter in enumerate(costaParameters):
	materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,component, parameter)

### Step 10: Create equation set ####################################################
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.ACTIVECONTRACTION, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

### Step 11: Create dependent field ##############################################
dependentField = CMISS.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

dependentField.CreateStart(dependentFieldUserNumber, region)
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Dependent")
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT)
dependentField.NumberOfVariablesSet(2)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN, 4)
for i in [1,2,3]:
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
# Initialise from undeformed geometry
for component in [1,2,3]:
    geometricField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,component, dependentField,CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)

dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -8.0)


### Step 12: Create independent field ############################################
independentField = CMISS.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber, independentField)
equationsSet.IndependentCreateFinish()



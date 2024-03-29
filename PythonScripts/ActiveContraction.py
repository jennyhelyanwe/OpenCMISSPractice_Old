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
	problemUserNumber,
    cellMLUserNumber) = range(1,18)

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
	materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
materialField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
materialField.CreateFinish()

# Set costa material parameters. 
MRParameters = [2.0, 6.0]
for component, parameter in enumerate(MRParameters):
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
dependentField.NumberOfVariablesSet(4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U1, 6) # storing strain
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U2, 6) # storing stress

for component in range (1,7):
    dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1, component, CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2, component, CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)

for i in range (1,5):
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U1, i, 1)
for i in range (1,7):
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U2, i, 1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i, 1)    
dependentField.CreateFinish()

# Initialise from undeformed geometry
for component in [1,2,3]:
    geometricField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,component, dependentField,CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component)

dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -8.0)

### Step 12: Create CellML environment ###########################################
ActiveContractionModelIndex = 1
cellML = CMISS.CellML()
cellML.CreateStart(cellMLUserNumber, region)
cellML.ModelImport("active_contraction.cellml", ActiveContractionModelIndex)

strain = ["E11", "E12", "E13", "E22", "E23", "E33"]
stressCauchy = ["Tdev11", "Tdev12", "Tdev13", "Tdev22", "Tdev23", "Tdev33"]

for  i in range (0,6):
    cellML.VariableSetAsKnown(ActiveContractionModelIndex, "equations"+strain[i])



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

### Step 8: Fibre field ###############################################################
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre")
fibreField.NumberOfVariablesSet(1)
fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)
fibreField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
for component in range (1,4):
	fibreField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, component, CMISS.FieldInterpolationTypes.CONSTANT)
fibreField.CreateFinish()
# Initialise the fibre rotation angles in radians
fibreFieldAngle = 30 * pi /180  # 30 degrees in radians
fibreAngles = [fibreFieldAngle, 0, 0]
for component, fibreAngle in enumerate(fibreAngles, 1):
	fibreField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, component, fibreAngle)

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
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 0.2)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2, 30.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 3, 12.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, 14.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 5, 14.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 5, 10.0)
materialField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 5, 18.0)




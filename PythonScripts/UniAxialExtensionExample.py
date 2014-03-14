######### This script attempts to solve a uniaxial extension, tri-linear lagrange, isotropic cube ###########################
## JZW 05/03/2014

import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings', 'python')))

### Initialise OpenCMISS
from opencmiss import CMISS

### Set problem parameters
height = 1.0
width = 1.0
length = 1.0

numberGlobalXElements = 1
numberGlobalYElements = 1
numberGlobalZElements = 1
numOfXi = 3

### Set arbitrary usernumbers which are unique to each object. 
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber,
    fibreFieldUserNumber,
    materialFieldUserNumber) = range(1,15)


### Set up diagnostics/debug
CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.IN, [1,2,3,4,5], "Diagnostics", ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

### Get computational node information for parallel computing
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet() # Get total rank
computationalNodeNumber = CMISS.ComputationalNodeNumberGet() # Get self rank number

### Set up 3D RC co-ordinate system. 
coordinateSystem = CMISS.CoordinateSystem()	# Initialise instance of CS
coordinateSystem.CreateStart(coordinateSystemUserNumber)	# Create CS
coordinateSystem.dimension = 3;	# 3D CS set
coordinateSystem.CreateFinish()	

### Create a region.
region = CMISS.Region()	# Initialise an instance of region which encapsulates simulation
region.CreateStart(regionUserNumber, CMISS.WorldRegion)	# The region is define as a daughter region of WorldRegion
region.label = "Region"	# The label of the region is named after the problem being solved. 
region.coordinateSystem = coordinateSystem	# Assign a coordinate system to the region. 
region.CreateFinish()

### Create a tri-linear lagrange basis
basis = CMISS.Basis()
basis.CreateStart(basisUserNumber)	# Create an instance of basis object identified by a user number. 
basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP	# Set its type to Lagrange_Hermite
basis.numberOfXi = numOfXi	# Set the dimensions of element co-ordinate system 
basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numOfXi #Set interpolation type in each xi direction. 
basis.quadratureNumberOfGaussXi = [2]*numOfXi	# Set the number of gauss points per xi
basis.CreateFinish()

### Create a mesh
mesh = CMISS.Mesh() # Initialise an instance of mesh
mesh.CreateStart(meshUserNumber, region, numOfXi) # Create mesh in the specified region with 3 elemental directions. 
mesh.NumberOfComponentsSet(1)	# Set the number of components of the mesh
mesh.NumberOfElementsSet(1)	# Set the total number of elements in the mesh

nodes = CMISS.Nodes()  # Initialise an instance of nodes object
nodes.CreateStart(region, 8)	# Create 8 nodes
nodes.CreateFinish() # Finish creation

elem = CMISS.MeshElements() # Initialise an instance of mesh elements. 
elem.CreateStart(mesh, 1, basis)	# Create an element in mesh component 1. 
elem.NodesSet(1,[1,2,3,4,5,6,7,8])	# The nodes in global element 1 are given user node numbers 1 to 8. 
elem.CreateFinish()

mesh.CreateFinish()

### Create a decomposition for the mesh 
decomposition = CMISS.Decomposition() # Initialise
decomposition.CreateStart(decompositionUserNumber, mesh) # Create decomposition for the mesh previously defined. 
decomposition.type = CMISS.DecompositionTypes.CALCULATED # Set the type of decomposition
decomposition.NumberOfDomainsSet(numberOfComputationalNodes) # There are as many domains as there are processors. 
decomposition.CreateFinish()

### Create a field for geometry
geometricField = CMISS.Field() # Initialise
geometricField.CreateStart(geometricFieldUserNumber, region) # Creates field in the defined region. 
geometricField.TypeSet(CMISS.FieldTypes.GEOMETRIC) # Set the field as geometric type. 
geometricField.MeshDecompositionSet(decomposition) # Set the mesh decomposition of field. 
geometricField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Geometry") # Set type to u, not the partial derivative or anything else. 
# Set the field components all to the single mesh component: i.e. the geometry is interpolated in the same fashion in all directions. 
fieldVarType = CMISS.FieldVariableTypes.U
geometricField.ComponentMeshComponentSet(fieldVarType,1,1) # Set the field component 1 for mesh component number 1
geometricField.ComponentMeshComponentSet(fieldVarType,2,1) # Set the field component 2 for mesh component number 1
geometricField.ComponentMeshComponentSet(fieldVarType,3,1) # Set the field component 3 for mesh component number 1
geometricField.CreateFinish()

### Update the geometric field parameters manually
fieldParSetType = CMISS.FieldParameterSetTypes.VALUES
geometricField.ParameterSetUpdateStart(fieldVarType, fieldParSetType)
# Node 1
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,1,1,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,1,2,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,1,3,0.0)
# Node 2
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,2,1,length)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,2,2,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,2,3,0.0)
# Node 3
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,3,1,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,3,2,width)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,3,3,0.0)
# Node 4
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,4,1,length)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,4,2,width)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,4,3,0.0)
# Node 5
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,5,1,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,5,2,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,5,3,height)
# Node 6
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,6,1,length)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,6,2,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,6,3,height)
# Node 7
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,7,1,0.0)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,7,2,width)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,7,3,height)
# Node 8
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,8,1,length)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,8,2,width)
geometricField.ParameterSetUpdateNodeDP(fieldVarType, fieldParSetType, 1,1,8,3,height)

### Create a fibre field and attach it to the geometric field
fibreField = CMISS.Field() # Initialise
fibreField.CreateStart(fibreFieldUserNumber,region) # Create with unique fibre field user number. 
fibreField.TypeSet(CMISS.FieldTypes.FIBRE) # Set type of field to fibre
fibreField.MeshDecompositionSet(decomposition) # Set decomposition of the field.
fibreField.GeometricFieldSet(geometricField) # Apply the fibre field on the geometry just created. 
#fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 2)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre") # Set parameter type to U and label it Fibre. 
fibreField.CreateFinish()

### Create material field
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber,region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL) # Set type of field to material properties
materialField.MeshDecompositionSet(decomposition) # Set decomposition to field
materialField.GeometricFieldSet(geometricField) # Set geometry to field
materialField.VariableLabelSet(fieldParSetType, "Material") # Label field variable
materialField.NumberOfVariablesSet(1) #There is one variable set
materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U,2) # With 2 components
# Apply field to the mesh component
materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
materialField.CreateFinish()

### Set the stiffness coefficient. 
### Set the exponential constants c1 to c5 respectively
CMISS.Field.ComponentValuesInitialiseDP(materialField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 2.0) 
CMISS.Field.ComponentValuesInitialiseDP(materialField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2, 6.0) 

### Create the dependent field 
dependentField = CMISS.Field() # Initialise a dependent field
dependentField.CreateStart(dependentFieldUserNumber, region)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Dependent") # label field
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL) # Set type of field to Geometric General.
dependentField.MeshDecompositionSet(decomposition) # Set decomposition
dependentField.GeometricFieldSet(geometricField) # Set geometry
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT) # Dependent field
dependentField.NumberOfVariablesSet(2) # One set of variables, one for u and one for boundary flow. 
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN,4)
# Set field components to the single mesh component
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 1,1)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 2,1)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 3,1)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 1,1)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 2,1)
dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 3,1)
# Set the 4th component to being element based interpolation instead of node based, 
# this component is for the hydrostatic pressure. 
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, 4, CMISS.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.DELUDELN, 4, CMISS.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.CreateFinish()

### Initialise dependent field. 
# Copy geometric field components to that of the dependent field as initial undeformed 
# conformation. 
CMISS.Field.ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1)

CMISS.Field.ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2, dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 2)

CMISS.Field.ParametersToFieldParametersComponentCopy(
geometricField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 3, dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 3)

# Set hydrostatic pressure
CMISS.Field.ComponentValuesInitialiseDP(dependentField, CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4, -8.0)

### Create equations 
equationsSetField = CMISS.Field() # Equations are also in a field
equationsSet = CMISS.EquationsSet() # Initialise an equation set. 
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY, CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.MOONEY_RIVLIN,
equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Add the material field to the equation
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()

# Add the dependent field to the equation
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

# Create equations from the equation set
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE # Sparse matrix sufficient for FEM. 
equations.outputType = CMISS.EquationsOutputTypes.NONE # Set output of equations to none
equationsSet.EquationsCreateFinish()

### Define problem
problem = CMISS.Problem() # Initialise problem object
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY, CMISS.ProblemTypes.FINITE_ELASTICITY, CMISS.ProblemSubTypes.NONE) # Specify type of problem solved
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

### Create problem solver
nonLinearSolver = CMISS.Solver() 
linearSolver = CMISS.Solver() # Initialise a linear and a non-linear solver. 
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver) # initialise control loop for the problem, there is only one type of control loop identifier. 
nonLinearSolver.outputType = CMISS.SolverOutputTypes.PROGRESS # set output as progress
nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.FD) # calculate jacobian using finite differences. 
nonLinearSolver.NewtonLinearSolverGet(linearSolver) # set newton method linear solver. 
linearSolver.linearType = CMISS.LinearSolverTypes.DIRECT # use direct not iterative linear solver
problem.SolversCreateFinish()

### Create solver equations and add equations set to solver equations
solver = CMISS.Solver() # Initialise solver object
solverEquations = CMISS.SolverEquations() # Initialise solver equations object
problem.SolverEquationsCreateStart() # Initialise solver equations in problem
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations) # load equation set into problem solver equations
solverEquations.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()


### Prescribe BC using absolute nodal parameters
BC = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(BC)

# Set BC for extension of 10% in positive x direction. 
fixed = CMISS.BoundaryConditionsTypes.FIXED
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,1,1, fixed, 0.0) # Inputs are: field, variable type, version number, derivative number, node user number, component number, condition (free, fixed, neumann, dirichlet, pressure, etc.), value. 
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,3,1, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,5,1, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,7,1, fixed, 0.0)

BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,2,1, fixed, 0.1*length)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,4,1, fixed, 0.1*length)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,6,1, fixed, 0.1*length)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,8,1, fixed, 0.1*length)

# Set fixed zero displacement for front. 
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,1,2, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,2,2, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,5,2, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,6,2, fixed, 0.0)

# Set fixed zero displacement for bottom face. 
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,1,3, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,2,3, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,3,3, fixed, 0.0)
BC.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,1,4,3, fixed, 0.0)

# Finish
solverEquations.BoundaryConditionsCreateFinish()

### Solve problem
problem.Solve()

### Export results
fields = CMISS.Fields()
fields.CreateRegion(region)
fields.NodesExport("../Results/UniAxialExtension", "FORTRAN")
fields.ElementsExport("../Results/UniAxialExtension","FORTRAN")
fields.Finalise()




# Basic Rules:

## Program at interface
First describe the interface in a separate class, then implement.

## Use templates when
- ...

## Use inheritance when
- ...

## Naming conventions
* Variable:
  - lower case with underscores
  - single letter when possible
  - if first letter, put that word in comment: n // name
  - if another letter, put that word with letter in [...]: a // n[a]me
  - up to 4 letters if possible
  - private variables end with underscore: `a_`;
    exceptions: m (Mesh), var (Vars)
* Template argument:
  - CamelCase for types and lower case for values
  - if alias needed, ends with underscore: `Mesh_`
  - else one letter
* Filename:
  - lower case, extensions h and cpp
* Class
  - ...
* No names starting with underscore, reserved
* For indices:
  `i` - generic 
  `c` - IdxCell
  `f` - IdxFace 
  `w` - MIdx
  `d` - direction
  `d` - index 0..dim
* Other:
  `u` - generic field
  `q` - neighbour cell id
  `o` - other
* Treat pointers same as values

## Other

* use Vect(0) and MIdx(0) instead of Vect::kZero

## Formatting
* initialization list: 4 spaces indent `:` and `,` starting each line.

## Toolbox classes
* combines multiple templated functions to avoid repeating template parameters
* acts as namespace
* contains only static functions

## Dereference pointers to references or values if possible
But: use pointer arguments to prevent passing an rvalue

## Lower bound for template argument:
  template <class B /*: A*/>
means that argument B needs to be a subtype of A
 
## Comments
- capitalized descriptive for a class or function 
  (e.g. Creates instance)
- capitalized imperative for expressions in implementation 
  (e.g. Create instance)
- non-capitalized for declarations 
  (e.g. buffer index)

## Data from Kernel returned by reference if possible

## Units

Example: `Scal s; // source [density/time]`

* length
* time
* mass
* velocity (length/time)
* area (length^2)
* volume (length^3)
* density (mass/volume)
* value (user defined, e.g. 1 for concentration)

## Member classes

* A class taking mesh as template argument 
should call it `M_` and define member class `M`.
* A class taking kernel as template argument
should call it `K_` and define member class `K`.
* A derived class can optionally define member class `P` for parent.
* A mesh class should define member classes
`Scal` and `Vect` and member variable `dim`.

### Ordering

* parent `P`
* kernel `K`
* mesh `M`
* constructor parameters `Par`
* dimension `dim`

### Dependencies and interfaces

Instances: Mesh=MeshStructured, Kernel=Hydro, Distr=Cubism

Interfaces:
- Mesh: cell connectivity, geometry.
  Must be a template argument for performance.
  Aware of: FieldCell, IdxCell
- Kernel: Run(), Comm(), Reduce(), ReadBuffer(), WriteBuffer(), GetComm()
  Aware of: Mesh (template), FieldCell, IdxCell, BlockInfo
- KernelFactory: Creates Kernel
  Aware of: Kernel, Mesh (template)
- Distr: Takes KernelFactory, instantiates Kernel for every block,
  reads GetComm(), implements communication
 
Requirements:
- Cubism should not depend on Hydro (only on Kernel and KernelFactory)
- Hydro should not depend on Cubism (only on Distr)
- Hydro is Kernel, HydroFactory is KernelFactory
- Ideally, Hydro doesn't depend on implementation of Mesh
- Cubism is Distr
- HydroFactory creates Hydro for a block described in Cubism-independent way
- Hydro accepts Mesh initalized by HydroFactory 
  and is not aware of Distr,
  but it (as well as other mesh-related functions like Gradient())
  aware of Mesh which has some distributed computing primitives 
  like Comm() and Reduce() (also Solve() to solve linear systems
  but that should be avoided)
- Distributed primivites in Mesh are not dependent on Cubism or even Distr.
  Typically, it is just a list of fields for communication.
- Why doesn't Mesh depend on Distr?
  This way Mesh is simpler and has less dependencies.
  Possible drawbacks: information from MPI is not available.
  All routines (e.g. solve linear system) rely on rigid interface
  and need to be performed from the outside
  (though this is necessary for any blocking setup).
- Interface of Hydro: 
    vector<FieldCell<Scal>*> GetComm()
    M& GetMesh()
- Interface of Cubism:
    void Step()
  Cubism (and any Distr) acts like a framework.
  But every Distr has a way to access and communicate data in blocks
  (in Cubism, it is BlockLab, Block_t and BlockInfo)
- The question is how to organize interaction between Distr and Kernel
  (in particular, Cubism and Hydro). Options:
  1) sufficient interface of Kernel
  How: same WriteBuffer() and ReadBuffer() but with some generic block buffer.
  Cons: too rigid and may require data copying
  2) some entity aware of both Hydro (implementation of Kernel)
  and Cubsm (implementation of Distr), different for every pair.
  How: visitor?
  3) make Cubism aware of Hydro or at least Kernel<MeshStructured>


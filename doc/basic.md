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
    exceptions: m (Mesh), par (Vars)
* Template argument:
  - CamelCase for types and lower case for values
  - if alias needed, ends with underscore: `Mesh_`
  - else one letter
* Filename:
  - lower case, extensions h and cpp
* Class
  - ...
* No names starting with underscore, reserved

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




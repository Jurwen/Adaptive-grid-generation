# Adaptive Grid Generation

The `isosurfacing` tool performs Longest Edge Bisection Refinement to generate an adaptive grid from an initial mesh and an implicit function, and by incorporating a [robust isosurfacing method](https://github.com/duxingyi-charles/Robust-Implicit-Surface-Networks/tree/main), we extract the surfaces for implicit complexes.

## Build

Use the following command to build: 

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
The program `isosurfacing` will be generated in the build file. 


## Usage

To use the `isosurfacing` tool, you must provide an initial mesh file and implicit function file as arguments, along with any desired options.

```bash
./isosurfacing <mesh> <function> [OPTIONS]
```

### Positional Arguments

- `mesh` : The path to the initial mesh file that will be used for isosurfacing. This file can either be a `.msh` or `.json` file. 
Examples of mesh files can be found in the `data/mesh` directory.
- `function` : The path to the implicit function file that is used to evaluate the isosurface.

### Options

- `-h, --help` : Show the help message and exit the program.
- `-t, --threshold` : Set the threshold value for the isosurface generation. This is a `DOUBLE` value that defines the precision level of the isosurfacing.
- `-o, --option` : Set the type of the implicit complexes from Implicit Arrangement(IA), Material Interface(MI), or Constructive Solid Geometry(CSG). This is a `STRING` value that takes input from "IA", "MI", and "CSG". The default type is "IA".
- `--tree` : The path to the CSG tree file that defines the set of boolean operations on the function. 
- `-c, --curve_network` : Set the switch of extracting only the Curve Network. Notice that Curve Network of all the above implicit complexes are different given the same set of input functions. This is a `BOOLEAN` type that takes in 1 or 0.
- `-m, --max-elements` : Set the maximum number of elements in the mesh after refinement. This is an `INT` value that limits the size of the generated mesh. If this value is a **negative** number, the mesh will be refined until the threshold value is reached.
- `-s,--shortest-edge` : Set the shortest length of edges in the mesh after refinement. This is a `DOUBLE` value that defines the shortest edge length.

## Example

The following is an example of how to use the `isosurfacing` tool with all available options:

```bash
./isosurfacing data/mesh/cube6.msh data/function/sphere.json -t 0.01 -o "IA" -m 10000 -s 0.05
```

In this example, `cube6.msh` is the initial mesh file, `sphere.json` is the implicit function file, `0.01` is the threshold value, "IA" is the type of implicit complexes, `10000` is the maximum number of elements, and `0.05` is the shortest edge length.

## Help

You can always run `./isosurfacing -h` to display the help message which provides usage information and describes all the options available.

## After Grid Generation

We have an off-the-shelf algorithm that extracts the isosurfacs robustly from the grid for implicit complexes.  

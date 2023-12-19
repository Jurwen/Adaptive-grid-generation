# Isosurfacing Tool

The `isosurfacing` tool performs Longest Edge Bisection Refinement to generate an adaptive grid from an initial mesh and an implicit function.

## Usage

To use the `isosurfacing` tool, you must provide an initial mesh file and implicit function file as arguments, along with any desired options.

```bash
./isosurfacing <mesh> <function> [OPTIONS]
```

### Positional Arguments

- `mesh` : The path to the initial mesh file that will be used for isosurfacing.
- `function` : The path to the implicit function file that is used to evaluate the isosurface.

### Options

- `-h, --help` : Show the help message and exit the program.
- `-t, --threshold` : Set the threshold value for the isosurface generation. This is a `FLOAT` value that defines the precision level of the isosurfacing.
- `-m, --max-elements` : Set the maximum number of elements in the mesh after refinement. This is an `INT` value that limits the size of the generated mesh. If this value is a negative integer, the mesh will be refined until the threshold value is reached.

## Example

The following is an example of how to use the `isosurfacing` tool with all available options:

```bash
./isosurfacing data/mesh/cube6.msh data/function/sphere.json -t 0.01 -m 10000 
```

In this example, `cube6.msh` is the initial mesh file, `sphere.json` is the implicit function file, `0.01` is the threshold value, and `10000` is the maximum number of elements.

## Help

You can always run `./isosurfacing -h` to display the help message which provides usage information and describes all the options available.

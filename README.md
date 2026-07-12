# Project2

A C++ project implementing **approximate/exact nearest neighbor search** and **clustering** algorithms, including:

- **LSH (Locality Sensitive Hashing)** for approximate nearest neighbor and range search
- **Hypercube projection** for approximate nearest neighbor and range search
- **Clustering** using various initialization, assignment, and update methods
- Support for **discrete Frechet / continuous curve** operations (`Curve`)

## Overview

The repository provides multiple executables that operate on a dataset of vectors or curves, allowing you to:

- Search for nearest neighbors using LSH or Hypercube methods
- Perform clustering on a dataset and evaluate the results (e.g., via silhouette score)
- Run exact vs. approximate search comparisons

## Project Structure

```
.
├── header/          # Header files (.h) for all data structures and algorithms
│   ├── Curve.h
│   ├── Dataset.h
│   ├── HashTable.h
│   ├── Hypercube.h
│   ├── InputParser.h
│   ├── Point.h
│   ├── clustering.h
│   ├── neighbour_search.h
│   └── vector_operations.h
├── src/              # Implementation files (.cpp)
│   ├── Curve.cpp
│   ├── Dataset.cpp
│   ├── HashTable.cpp
│   ├── Hypercube.cpp
│   ├── clustering.cpp
│   ├── neighbour_search.cpp
│   └── vector_operations.cpp
├── src_main/         # Executable entry points
│   ├── cluster_main.cpp   # Clustering executable
│   ├── hcube_main.cpp     # Hypercube-based search executable
│   ├── lsh_main.cpp       # LSH-based search executable
│   └── search_main.cpp    # General search executable
├── unit_tests/       # Unit tests (using the acutest framework)
│   ├── test.cpp
│   └── utilities.h
└── FredFilesNeeded/  # Additional supporting files
```

## Building

The project is written primarily in **C++** (with a small amount of **C**). It can be compiled with a standard C++ compiler such as `g++`.

Example (adjust flags/paths as needed for your setup):

```bash
g++ -std=c++17 -Iheader src/*.cpp src_main/cluster_main.cpp -o cluster
g++ -std=c++17 -Iheader src/*.cpp src_main/lsh_main.cpp -o lsh
g++ -std=c++17 -Iheader src/*.cpp src_main/hcube_main.cpp -o hcube
g++ -std=c++17 -Iheader src/*.cpp src_main/search_main.cpp -o search
```

> A `Makefile` is recommended if not already present, to simplify building all executables at once.

## Usage

Each executable typically accepts input parameters via the command line, such as:

- `-i` : input file (dataset)
- `-q` : query file
- `-o` : output file
- `-k`, `-L` : LSH parameters (number of hash functions / hash tables)
- `-M`, `-probes` : Hypercube parameters

Refer to each `*_main.cpp` file in `src_main/` for the exact set of supported arguments.

### Examples

```bash
./lsh -i input.dat -q query.dat -o output.txt -k 4 -L 5
./hcube -i input.dat -q query.dat -o output.txt -M 10 -probes 2
./cluster -i input.dat -c cluster_config.txt -o output.txt
```

## Testing

Unit tests are located in `unit_tests/` and use the lightweight [acutest](https://github.com/mity/acutest) testing framework.

```bash
g++ -std=c++17 -Iheader -Iunit_tests src/*.cpp unit_tests/test.cpp -o run_tests
./run_tests
```

## Language Composition

- **C++**: 99.2%
- **C**: 0.8%

## License

No license file is currently included in this repository. Please contact the repository owner for usage terms.

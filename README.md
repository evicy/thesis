## Identifying Clusters in Graph Representations of Genomes
Dynamic programming algorithm to find maximum-scoring disjoint paths in n-layered bubble graphs.

### Compilation
```
cd main_program/build
make main
```

### Usage
The program accepts Elastic-degenerate strings, EDS for short, which represents pan-genomes. To generate a EDS input, you can use [EDSO](https://github.com/webmasterar/edso).

```
./main generated_eds_string
```

### Tests
The tests use GoogleTests, make sure you have the module installed and add the path to CMakeLists.txt file.
```
make test
./test
```
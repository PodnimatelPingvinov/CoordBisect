### coordbisect.cpp
Takes 4 command line arguments:
```sh
mpirun -np <proc_num> coordbisect <k> <n1> <n2> <result.bin>
```
Where:
- `k` — number of domains.
- `n1` — lines number in mesh.
- `n2` — columns number in mesh.
- `result.bin` — resulting file, which will contain decomposed mesh in binary form.

### view.cpp
Program for viewing resulting binary file in human-readable form. Usage:
```sh
./view <result.bin>
```
It will print in standart output *result.bin* content. Each line describes one mesh node in
the following foramt:
```sh
i j x y domain
```
Where:
- `i` — line number.
- `j` — column number.
- `x` — x coordinate.
- `y` — y coordinate.
- `domain` — domain number of this node.

# Technology-Dependent Synthesis and Optimization of Circuits for Small S-boxes

This is a tool designed for the hardware implementation of small S-boxes. It utilizes the Breadth First Search algorithm and a heuristic pruning strategy. For more details, please refer to the paper \[1\](https://cic.iacr.org/p/1/4/20) and \[2\](https://eprint.iacr.org/2025/103).

## Functionalities

Given an *n*-bit S-box with some specific configuration files and parameters, this tool can generate a circuit with small area.

The S-box supports various sizes, including 3-bit, 4-bit, 5-bit, and 6-bit S-boxes. However, generally speaking, this tool is most effective for all 3-bit S-boxes, most 4-bit S-boxes, and a few 5-bit S-boxes.

The S-box does not need to be bijective.


## Usage

There are four versions of the code available:
- `full_single.cpp`: A single-threaded version without pruning strategy.
- `full_multi.cpp`: A multi-threaded version without pruning strategy.
- `reduced_single.cpp`: A single-threaded version with pruning strategy.
- `reduced_multi.cpp`: A multi-threaded version with pruning strategy.

After several experiments we find that the multi-threaded versions can hardly improve the computational. For small S-boxes, it is recommended to use `full_single.cpp` or `full_multi.cpp`. For large S-boxes, `reduced_single.cpp` or `reduced_multi.cpp` is a better choice.

The usage steps are as follows:

### Step 1: Configure the S-box size
Modify the 26th line of the code according to the S-box size:
- For a 3-bit S-box, use `#define UINT uint8_t`.
- For a 4-bit S-box, use `#define UINT uint16_t`.
- For a 5-bit S-box, use `#define UINT uint32_t`.
- For a 6-bit S-box, use `#define UINT uint64_t`.

### Step 2: Compile the code
Use the `g++` compiler in the system terminal to generate an executable file. Taking the `reduced_single.cpp` version as an example:
- On Windows:
```bash
g++ -fopenmp -Wall -g -O2 -static-libgcc -std=c++17 reduced_single.cpp -o reduced_single.exe
```
- On Linux:
```bash
g++ -fopenmp -Wall -g -O2 -static-libgcc -std=c++17 reduced_single.cpp -o reduced_single
```

### Step 3: Run the executable file
Run the executable file in the system terminal to generate the circuit of the S-box. Taking the 3-way S-box as an example:
- On Windows:
```bash
.\reduced_single.exe --cipher "0702040501060300" --areaconf SMIC65.conf --depthconf= --number 8 --depth 8 --result "3-way_SMIC65_area.txt"
```
- On Linux:
```bash
./reduced_single --cipher "b4 c6 9a" --areaconf SMIC65.conf --depthconf= --number 8 --depth 8 --result "3-way_SMIC65_area.txt"
```
- `cipher`: Specifies the S-box. It can be either the LUT representation encoded in hexadecimal (e.g., `0702040501060300`) or the bitsliced representation (concatenation of the *n* value vectors of its coordinate functions, e.g., `b4 c6 9a`). This parameter is required.
- `areaconf`: Specifies the file containing the area information of the used gates. You can refer to the `SMIC65.conf` file or `SMIC130.conf` for the format (the name and the value should be separated by `=`, do not support whitespace, names are not case sensitive). This parameter is required. All available logic gates are shown as follows:
  - NOT
  - AND
  - NAND
  - NANDN
  - OR
  - NOR
  - NORN
  - XOR
  - XNOR
  - AND3
  - NAND3
  - NANDN3
  - OR3
  - NOR3
  - NORN3
  - XOR3
  - XNOR3
  - MUX
  - MUXI
  - AO21
  - AOI21
  - OA21
  - OAI21<br>
If a certain gate is not needed, delete the data of the corresponding line in the file. To learn more about gates, please refer to the paper.
- `depthconf`: Specifies the file containing the delay information of the used gates. This is an optional parameter. If not provided, the default delay information will be used, assuming the delay of all gates is 1. If you need to specify custom delay information, you can create a file and provide its name. You can refer to the `SMIC130_Depth.conf` file for the format (the name and the value should be separated by `=`, do not support whitespace, names are not case sensitive).
- `number`: Specifies the maximum number of used gates in the circuit implementation to be generated. The actual number of gates of the generated circuits will not exceed this value. A larger value may result in a more optimized result but will also increase the computational time.
- `depth`: Specifies the maximum delay of the generated circuit implementations. The actual delay of the generated circuits will not exceed this value.
- `result`: Specifies the output file name. The generated circuit information, including the area, delay, used gate types and used gate numbers, will be saved in this file. The input values of S-box are represented as `x0`, `x1`, `x2`, etc., and the output values of S-box are represented as `y0`, `y1`, `y2`, etc, where `x0` and `y0` represent the least significant bit.


## References

[1] Zihao Wei, Siwei Sun, Fengmei Liu, Lei Hu and Zhiyu Zhang.: Technology-Dependent Synthesis and Optimization of Circuits for Small S-boxes. https://cic.iacr.org/p/1/4/20.

[2] Zihao Wei, Siwei Sun, Fengmei Liu, Lei Hu and Zhiyu Zhang.: Technology-Dependent Synthesis and Optimization of Circuits for Small S-boxes. https://eprint.iacr.org/2025/103.

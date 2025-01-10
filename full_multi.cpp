#include <algorithm>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>

#include <time.h>

using namespace std;

/*****
 * if the cipher is a 3-bit S-box, define uint8_t as UINT
 * if it is a 4-bit S-box, define uint16_t as UINT
 * if it is a 5-bit S-box, use uint32_t
 * if it is a 6-bit S-box, use uint64_t
*****/
// define UINT type
#define UINT uint8_t

/*****
 * if a global variable name has the form like "gate_xx", it means the variable is one of the inherent properties of the gate
 * if a global variable name has the form like "opt_xx", it means the variable depends on the library file
/*****
 *  predefine all executable gates
 * if users want to use their own gate, add the gate name and its characteristics
*****/
// all possible gate names
vector<string> gate_name = {"NOT", "AND", "NAND", "NANDN", "OR", "NOR", "NORN", "XOR", "XNOR", "AND3", "NAND3", "NANDN3", "OR3", "NOR3", "NORN3", "XOR3", "XNOR3", "MUX", "MUXI", "AO21", "AOI21", "OA21", "OAI21"};
enum gate_value
{
    NotDefined,
    NOT,
    AND,
    NAND,
    NANDN,
    OR,
    NOR,
    NORN,
    XOR,
    XNOR,
    AND3,
    NAND3,
    NANDN3,
    OR3,
    NOR3,
    NORN3,
    XOR3,
    XNOR3,
    MUX,
    MUXI,
    AO21,
    AOI21,
    OA21,
    OAI21
};
unordered_map<string, gate_value> gate_map =
    {
        {"NOT", NOT},
        {"AND", AND},
        {"NAND", NAND},
        {"NANDN", NANDN},
        {"OR", OR},
        {"NOR", NOR},
        {"NORN", NORN},
        {"XOR", XOR},
        {"XNOR", XNOR},
        {"AND3", AND3},
        {"NAND3", NAND3},
        {"NANDN3", NANDN3},
        {"OR3", OR3},
        {"NOR3", NOR3},
        {"NORN3", NORN3},
        {"XOR3", XOR3},
        {"XNOR3", XNOR3},
        {"MUX", MUX},
        {"MUXI", MUXI},
        {"AO21", AO21},
        {"AOI21", AOI21},
        {"OA21", OA21},
        {"OAI21", OAI21}};
// the number of operators of gates
unordered_map<string, int> gate_number =
    {
        {"NOT", 1},
        {"AND", 2},
        {"NAND", 2},
        {"NANDN", 2},
        {"OR", 2},
        {"NOR", 2},
        {"NORN", 2},
        {"XOR", 2},
        {"XNOR", 2},
        {"AND3", 3},
        {"NAND3", 3},
        {"NANDN3", 3},
        {"OR3", 3},
        {"NOR3", 3},
        {"NORN3", 3},
        {"XOR3", 3},
        {"XNOR3", 3},
        {"MUX", 3},
        {"MUXI", 3},
        {"AO21", 3},
        {"AOI21", 3},
        {"OA21", 3},
        {"OAI21", 3}};
// the symmetry of operations, i.e. 1 means symmetric, 2 means asymmetric, 3 and 4 mean semi-symmetric, where NANDN3(a, b, c) = NANDN3(a, c, b) is the example of 3, and AO21(a, b, c) = AO21(b, a, c) is the example of 4
unordered_map<string, int> gate_symmetry =
    {
        {"NOT", 1},
        {"AND", 1},
        {"NAND", 1},
        {"NANDN", 2},
        {"OR", 1},
        {"NOR", 1},
        {"NORN", 2},
        {"XOR", 1},
        {"XNOR", 1},
        {"AND3", 1},
        {"NAND3", 1},
        {"NANDN3", 3},
        {"OR3", 1},
        {"NOR3", 1},
        {"NORN3", 3},
        {"XOR3", 1},
        {"XNOR3", 1},
        {"MUX", 2},
        {"MUXI", 2},
        {"AO21", 4},
        {"AOI21", 4},
        {"OA21", 4},
        {"OAI21", 4}};
// operations that the library file stands, which are the actual operations used in the program, the operations are sorted in ascending order based on their cost, 
vector<string> opt_name;
// cost of operations that the library file provides
unordered_map<string, int> opt_cost;
// depth of operations that the library file provides
unordered_map<string, int> opt_depth;
// number of operations used in the final circuit
map<string, int> opt_number;

// define a class
/*****
 * every circuit is a structure including 
 * the used gate amount, 
 * a vector recording the points in this circuit, 
 * a vector recording the depth corresponding to the points, 
 * the received target point amount, 
 * the creation points of the latest point
 * the creation operator
 * the previous circuit number
 * the number of this circuit in a vector
*****/
class circuit
{
public:
    // points information
    vector<UINT> point;
    vector<int> depth;
    int gate_number;
    int target_number;
    // path information
    vector<UINT> parent_point;
    string opt;
    int previous_Nu;
    int Nu;

    // initiate a circuit
    circuit() {}

    circuit(const vector<UINT> &P, const vector<int> &D) : point(P), depth(D)
    {
        gate_number = 0;
        target_number = 0;
        Nu = 0;
    }

    // generate a new circuit from the old one
    void generate(const circuit &old, UINT new_point, int new_depth, bool is_target, const vector<UINT> &parents, string opt_way)
    {
        // update point
        point.reserve(old.point.size() + 1);
        point = old.point;
        point.push_back(new_point);
        depth.reserve(old.depth.size() + 1);
        depth = old.depth;
        depth.push_back(new_depth);
        gate_number = old.gate_number + 1;
        target_number = old.target_number + is_target;
        // update path
        parent_point = parents;
        opt = opt_way;
        previous_Nu = old.Nu;
    }

    // extract the points whose depth less than or equal to allowable_depth and their corresponding depth
    void extract(int allowable_depth, vector<UINT> &result_point, vector<int> &result_depth)
    {
        result_point.reserve(point.size());
        result_depth.reserve(depth.size());
        for (vector<int>::size_type i = 0; i < depth.size(); ++i)
        {
            if (depth[i] <= allowable_depth)
            {
                result_point.push_back(point[i]);
                result_depth.push_back(depth[i]);
            }
        }
    }

    // judge if the new_point is repeated
    bool repeated_point(UINT new_point)
    {
        return find(point.begin(), point.end(), new_point) != point.end();
    }

    // judge if the circuit will never reach all targets with at most max_number gates
    bool unpromising_circuit(int max_number, int bit_number)
    {
        return bit_number - target_number + gate_number > max_number;
    }

    // record the information of how to generate the latest point
    string path(int bit_number, const vector<UINT> &Y)
    {
        if (opt_number.count(opt) == 0)
            opt_number[opt] = 1;
        else
            ++opt_number[opt];
        vector<string> parents;
        for (UINT ele : parent_point)
        {
            int index = find(point.begin(), point.end(), ele) - point.begin();
            stringstream ss;
            if (index < bit_number)
                ss << "x" << bit_number - 1 - index;
            else
                ss << "t" << index - bit_number + 1;
            parents.push_back(ss.str());
        }
        stringstream result;
        result << "t" << gate_number << " = " << opt << "(" << parents[0];
        for (vector<string>::size_type i = 1; i < parents.size(); ++i)
            result << ", " << parents[i];
        result << ")\t" << (double)(*depth.rbegin()) / 100;
        int index = find(Y.begin(), Y.end(), *point.rbegin()) - Y.begin();
        if (index < bit_number)
            result << "  [y" << bit_number - 1 - index << "]";
        return result.str();
    }

    // the collation of circuit
    bool operator<(const circuit &c) const
    {
        return (target_number > c.target_number) || ((target_number == c.target_number) && (gate_number < c.gate_number));
    }

    friend ostream &operator<<(ostream &output, const circuit &c)
    {
        output << "point: ";
        for (auto ele : c.point)
            output << (int)ele << " ";
        output << "\ndepth: ";
        for (auto ele : c.depth)
            output << ele << " ";
        output << "\ngate number: " << c.gate_number << endl;
        output << "target number: " << c.target_number << endl;
        output << "parents: ";
        for (auto ele : c.parent_point)
            output << (int)ele << " ";
        output << "\noperation: " << c.opt << endl;
        output << "from No." << c.previous_Nu << endl;
        output << "my No." << c.Nu << endl;
        return output;
    }
};

// record the aggregate of circuits according to the cost
unordered_map<int, vector<circuit>> cost_dic;
// record the aggregate of circuits according to the gate number, only have the information of point
unordered_map<int, set<vector<UINT>>> number_dic;

// parse the information of S-box
int parse_function(string str, vector<UINT> &X, vector<UINT> &Y)
{
    int bit_number;
    if (str.find(' ') != string::npos) // str is bit sliced representation
    {
        string temp;
        stringstream ss(str);
        bit_number = 0;
        while (getline(ss, temp, ' '))
        {
            Y.push_back((UINT)stoul(temp, nullptr, 16));
            ++bit_number;
        }
    }
    else // str is LUT representation
    {
        string temp;
        int l = str.length() / 2;
        bit_number = (int)(log(l) / log(2));
        vector<int> sbox(l);
        for (int i = 0; i < l; ++i)
        {
            temp = str.substr(2 * i, 2);
            sbox[i] = (UINT)stoul(temp, nullptr, 16);
        }
        for (int i = 0; i < bit_number; ++i)
        {
            UINT y = 0;
            for (int ind = 0; ind < l; ++ind)
                y = (y << 1) | ((sbox[ind] >> (bit_number - 1 - i)) & 1);
            Y.push_back(y);
        }
    }
    for (int i = 0; i < bit_number; ++i)
    {
        UINT x = 0;
        for (int ind = 0; ind < pow(2, bit_number); ++ind)
            x = (x << 1) | ((ind >> (bit_number - 1 - i)) & 1);
        X.push_back(x);
    }
    return bit_number;
}

// greatest common divisor
int gcd(int p, int q)
{
    int temp = p % q;
    while (temp != 0)
    {
        p = q;
        q = temp;
        temp = p % q;
    }
    return q;
}

// extract all different operators and their corresponding depth from configuration file, and return the greatest common divisor of gate cost
int extract_gate(string area_conf_name, string depth_conf_name)
{
    fstream conf_file;
    string line;
    multimap<int, string> opt_order;

    // extract information of cost
    conf_file.open(area_conf_name, ios::in);
    while (getline(conf_file, line))
    {
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        int len = line.find('=');
        string name = line.substr(0, len);
        if (find(gate_name.begin(), gate_name.end(), name) != gate_name.end())
        {
            int cost = (int)(stod(line.substr(len + 1, line.length() - len - 1)) * 100);
            opt_order.insert(pair<int, string>(cost, name));
            opt_cost[name] = cost;
            opt_depth[name] = 100;
        }
    }
    conf_file.close();

    // extract information of depth
    if (depth_conf_name != "")
    {
        conf_file.open(depth_conf_name, ios::in);
        while (getline(conf_file, line))
        {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            int len = line.find('=');
            string name = line.substr(0, len);
            if (opt_depth.find(name) != opt_depth.end())
            {
                int depth = (int)(stod(line.substr(len + 1, line.length() - len - 1)) * 100);
                opt_depth[name] = depth;
            }
        }
        conf_file.close();
    }

    // calculate the greatest common divisor
    opt_name.reserve(opt_order.size());
    int step = opt_order.begin()->first;
    for (auto ele : opt_order)
    {
        step = gcd(ele.first, step);
        opt_name.push_back(ele.second);
    }
    return step;
}

// calculate the number of combinations
// symmetry = 1: C_(total_number)^(choice_number)
// symmetry = 2: A_(total_number)^(choice_number)
// symmetry = 3 or 4: C_(total_number)^(choice_number) * (choice_number)
int number_of_combination(int total_number, int choice_number, int symmetry)
{
    int result = total_number;
    // symmetry = 2
    for (int i = total_number - 1; i > total_number - choice_number; --i)
        result *= i;
    if (symmetry == 1)
        for (int i = choice_number; i > 1; --i)
            result /= i;
    if (symmetry == 3 || symmetry == 4)
        for (int i = choice_number - 1; i > 1; --i)
            result /= i;
    return result;
}

// descending order
bool CMP(bool lhs, bool rhs)
{
    return lhs > rhs;
}

// push the result to result_point and result_depth
void display_result(const vector<UINT> &candidate_point, const vector<int> &candidate_depth, int symmetry, const vector<bool> &index, vector<vector<UINT>> &result_point, vector<int> &result_depth)
{
    vector<UINT> temp_point;
    int temp_depth = 0;
    for (vector<bool>::size_type i = 0; i < index.size(); ++i)
    {
        if (index[i])
        {
            temp_point.push_back(candidate_point[i]);
            if (candidate_depth[i] > temp_depth)
                temp_depth = candidate_depth[i];
        }
    }
    if (symmetry == 1)
    {
        result_point.push_back(temp_point);
        result_depth.push_back(temp_depth);
    }
    if (symmetry == 2)
    {
        sort(temp_point.begin(), temp_point.end());
        do
        {
            result_point.push_back(temp_point);
            result_depth.push_back(temp_depth);
        } while (next_permutation(temp_point.begin(), temp_point.end()));
    }
    if (symmetry == 3 || symmetry == 4)
    {
        vector<UINT>::size_type k;
        if (symmetry == 3)
            k = 0;
        else
            k = temp_point.size() - 1;
        sort(temp_point.begin(), temp_point.end());
        do
        {
            vector<vector<UINT>>::size_type i;
            for (i = 0; i < result_point.size(); ++i)
            {
                if (result_point[i][k] == temp_point[k])
                    break;
            }
            if (i == result_point.size())
            {
                result_point.push_back(temp_point);
                result_depth.push_back(temp_depth);
            }
        } while (next_permutation(temp_point.begin(), temp_point.end()));       
    }
}

// generate the combinations
void vector_of_combination(const vector<UINT> &candidate_point, const vector<int> &candidate_depth, int choice_number, int symmetry, vector<vector<UINT>> &result_point, vector<int> &result_depth)
{
    int total_number = (int)candidate_point.size();
    int number = number_of_combination(total_number, choice_number, symmetry);
    result_point.reserve(number);
    result_depth.reserve(number);
    vector<bool> index(total_number, false);

    for (int i = 0; i < choice_number; ++i)
        index[i] = true;
    display_result(candidate_point, candidate_depth, symmetry, index, result_point, result_depth);
    for (int i = 0; i < total_number - 1; ++i)
        if (index[i] && !index[i + 1])
        {
            index[i] = false;
            index[i + 1] = true;
            sort(index.begin(), index.begin() + i, CMP);
            display_result(candidate_point, candidate_depth, symmetry, index, result_point, result_depth);
            i = -1;
        }
}

/*****
 *  predefine all executable gates
 * if users want to use their own gate, add the expression
*****/
// generate new point according to parent points and gate
UINT operation(string opt, const vector<UINT> &parent)
{
    switch (gate_map[opt])
    {
    case NOT:
        return ~parent[0];
    case AND:
        return (parent[0] & parent[1]);
    case NAND:
        return ~(parent[0] & parent[1]);
    case NANDN:
        return ~(~parent[0] & parent[1]);
    case OR:
        return (parent[0] | parent[1]);
    case NOR:
        return ~(parent[0] | parent[1]);
    case NORN:
        return ~(~parent[0] | parent[1]);
    case XOR:
        return (parent[0] ^ parent[1]);
    case XNOR:
        return ~(parent[0] ^ parent[1]);
    case AND3:
        return (parent[0] & parent[1] & parent[2]);
    case NAND3:
        return ~(parent[0] & parent[1] & parent[2]);
    case NANDN3:
        return ~(~parent[0] & parent[1] & parent[2]);
    case OR3:
        return (parent[0] | parent[1] | parent[2]);
    case NOR3:
        return ~(parent[0] | parent[1] | parent[2]);
    case NORN3:
        return ~(~parent[0] | parent[1] | parent[2]);
    case XOR3:
        return (parent[0] ^ parent[1] ^ parent[2]);
    case XNOR3:
        return ~(parent[0] ^ parent[1] ^ parent[2]);
    case MUX:
        return ((parent[0] & parent[1]) ^ (~parent[0] & parent[2]));
    case MUXI:
        return ~((parent[0] & parent[1]) ^ (~parent[0] & parent[2]));
    case AO21:
        return ((parent[0] & parent[1]) | parent[2]);
    case AOI21:
        return ~((parent[0] & parent[1]) | parent[2]);
    case OA21:
        return ((parent[0] | parent[1]) & parent[2]);
    case OAI21:
        return ~((parent[0] | parent[1]) & parent[2]);
    default:
        return parent[0];
    }
}

// judge if the vector of point is repeated, if it is, return true, else update the number_dic and return false
bool repeated_vector_point(int number_key, const vector<UINT> &new_vector)
{
    bool result;
    #pragma omp critical 
    {
        auto iter = number_dic.find(number_key);
        if (iter != number_dic.end() && iter->second.find(new_vector) != iter->second.end())
            result = true;
        else
        {
            if (iter == number_dic.end())
            {
                set<vector<UINT>> number_value = {new_vector};
                number_dic[number_key] = number_value;
            }
            else
                iter->second.insert(new_vector);
            result = false;
        }
    }
    return result;
}

// if all target points are achieved, return total cost, if not, return 0
int Search(int max_number, int max_depth, int step, const vector<UINT> &Y, int bit_number)
{
    int cost_key = step;
    int result = 0;
    int min_depth = opt_depth[opt_name[0]];

    for (auto ele : opt_depth)
        if (ele.second < min_depth)
            min_depth = ele.second;
    
    while (cost_key <= max_number * opt_cost[*(--opt_name.end())])
    {
        // generate a list contains all circuits whose cost is cost_key
        list<circuit> cost_value_list;
        for (string Gate : opt_name)
        {
            auto cost_dic_iter = cost_dic.find(cost_key - opt_cost[Gate]);
            if (cost_dic_iter != cost_dic.end())
                for (circuit element : cost_dic_iter->second)
                {
                    vector<UINT> valid_point;
                    vector<int> valid_depth;
                    element.extract(max_depth - opt_depth[Gate], valid_point, valid_depth);
                    vector<vector<UINT>> temp_point;
                    vector<int> temp_depth;
                    vector_of_combination(valid_point, valid_depth, gate_number[Gate], gate_symmetry[Gate], temp_point, temp_depth);
                    #pragma omp parallel for
                    for (vector<UINT>::size_type i = 0; i < temp_point.size(); ++i)
                    {
                        #pragma omp flush(result)
                        if (!result)
                        {
                            int new_point = operation(Gate, temp_point[i]);
                            int new_depth = temp_depth[i] + opt_depth[Gate];
                            bool is_target = (find(Y.begin(), Y.end(), new_point) != Y.end());
                            // judge if the point is valid
                            if (!(element.repeated_point(new_point) || ((new_depth == max_depth && !is_target) || (new_depth < max_depth && !is_target && new_depth + min_depth > max_depth))))
                            {
                                circuit new_circuit;
                                new_circuit.generate(element, new_point, new_depth, is_target, temp_point[i], Gate);
                                vector<UINT> sorted_point(new_circuit.point.begin() + bit_number, new_circuit.point.end());
                                sort(sorted_point.begin(), sorted_point.end());
                                // judge if the circuit is valid
                                if (!(new_circuit.unpromising_circuit(max_number, bit_number) || repeated_vector_point(new_circuit.gate_number, sorted_point)))
                                {
                                    if (new_circuit.target_number == bit_number)
                                    {
                                        vector<circuit> cost_value_vector = {new_circuit};
                                        cost_dic[cost_key] = cost_value_vector;
                                        result = cost_key;
                                        #pragma omp flush(result)
                                    }
                                    else
                                    {
                                        #pragma omp critical 
                                        cost_value_list.push_back(new_circuit);
                                    }
                                }
                            }
                        }
                    }
                    if (result)
                        return result;
                }
        }

        // transform the list to the ordered vector, and store the vector in the cost_dic
        if (!cost_value_list.empty())
        {
            cost_value_list.sort();
            vector<circuit> cost_value_vector(cost_value_list.begin(), cost_value_list.end());
            for (vector<circuit>::size_type i = 0; i < cost_value_vector.size(); ++i)
                cost_value_vector[i].Nu = i;
            cost_dic[cost_key] = cost_value_vector;
        }

        // update cost_key
        cost_key += step;
    }
    return result;
}

int main(int argc, char *argv[])
{
    int opt;
    const char *shortopts = "";
    static struct option longopts[] =
        {
            {"cipher", required_argument, NULL, 'c'},
            {"areaconf", required_argument, NULL, 'f'},
            {"depthconf", optional_argument, NULL, 'F'},
            {"number", required_argument, NULL, 'n'},
            {"depth", required_argument, NULL, 'd'},
            {"result", required_argument, NULL, 'r'},
            {0, 0, 0, 0}};
    int bit_number, max_number;
    double max_depth;
    vector<UINT> X, Y;
    string area_conf_name, depth_conf_name, result_name, name;
    fstream result_file;

    while ((opt = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1)
    {
        switch (opt)
        {
        case 'c': // get the bit number and the bit-sliced form of input and output of S-Box according to "cipher"
            bit_number = parse_function(optarg, X, Y);
            break;
        case 'f': // get the area configure file name according to "areaconf"
            area_conf_name = optarg;
            break;
        case 'F': // get the depth configure file name according to "depthconf"
            depth_conf_name = optarg;
            if (depth_conf_name == "")
                name = "default configure";
            else
                name = depth_conf_name;
            break;
        case 'n': // get the maximum number of gate of circuit according to "number"
            max_number = stoi(optarg);
            break;
        case 'd': // get the maximun depth of circuit according to "depth"
            max_depth = stod(optarg);
            break;
        case 'r': // get the result file name according to "result"
            result_name = optarg;
            break;
        default:
            cout << "Invalid parameter" << endl;
            break;
        }
    }

    if (pow(2, bit_number) / 8 != sizeof(UINT))
    {
        cout << "Error: data type does not match the requirement of S-box, please redefine UINT" << endl;
        return 0;
    }

    result_file.open(result_name, ios::out);
    result_file << "Start points:" << endl;
    for (int i = 0; i < bit_number; ++i)
        result_file << (int)X[i] << " ";
    result_file << "\nTarget points:" << endl;
    for (int i = 0; i < bit_number; ++i)
        result_file << (int)Y[i] << " ";
    result_file << "\nArea configure file: " << area_conf_name << ", depth configure file: " << name << "\nMaximum number: " << max_number << ", maximum depth: " << max_depth << endl << endl;

    // get all possible gates according to area_conf_name, and get the step
    int step = extract_gate(area_conf_name, depth_conf_name);

    circuit start(X, vector<int>(bit_number, 0));
    vector<circuit> V = {start};
    cost_dic[0] = V;

    clock_t t1 = clock();
    int cost = Search(max_number, (int)(max_depth * 100), step, Y, bit_number);
    clock_t t2 = clock();
    cout << (double)(t2 - t1) / 1000 << endl;

    if (cost)
    {
        int Max_Depth = 0;
        for (int ele : cost_dic[cost][0].depth)
        {
            if (ele > Max_Depth)
                Max_Depth = ele;
        }
        cout << "success\nTotal cost: " << (double)cost / 100 << "\nMax depth: " << (double)Max_Depth / 100 << endl;
        result_file << "success\nTotal cost: " << (double)cost / 100 << "\nMax depth: " << (double)Max_Depth / 100 << endl;
        vector<string> path;
        int Nu = 0;
        while (cost)
        {
            circuit process = cost_dic[cost][Nu];
            path.push_back(process.path(bit_number, Y));
            cost -= opt_cost[process.opt];
            Nu = process.previous_Nu;
        }
        reverse(path.begin(), path.end());
        for (string ele : path)
        {
            cout << ele << endl;
            result_file << ele << endl;
        }
        result_file << endl;
        for (auto ele : opt_number)
            result_file << ele.first << ": " << ele.second << endl;
    }
    else
    {
        cout << "fail" << endl;
        result_file << "fail" << endl;
    }

    /*
    circuit start(X, vector<int>(bit_number, 0));
    vector<UINT> P;
    P.push_back(1);
    P.push_back(2);
    string ss = "AND";
    // first situation
    list<circuit> L;
    clock_t s1 = clock();
    for (UINT i = 0; i < 0xfff; ++i)
    {
        circuit second;
        bool b = find(Y.begin(), Y.end(), i) != Y.end();
        second.generate(start, i, 2, b, P, ss);
        L.push_back(second);
    }
    vector<circuit> V1(L.begin(), L.end());
    sort(V1.begin(), V1.end());
    clock_t e1 = clock();
    cout << (e1 - s1) << endl;

    // second situation
    multiset<circuit> S;
    clock_t s2 = clock();
    for (UINT i = 0; i < 0xfff; ++i)
    {
        circuit second;
        bool b = find(Y.begin(), Y.end(), i) != Y.end();
        second.generate(start, i, 2, b, P, ss);
        S.insert(second);
    }
    vector<circuit> V2(S.begin(), S.end());
    clock_t e2 = clock();
    cout << (e2 - s2) << endl;

    // third situation
    list<circuit> newL;
    clock_t s3 = clock();
    for (UINT i = 0; i < 0xfff; ++i)
    {
        circuit second;
        bool b = find(Y.begin(), Y.end(), i) != Y.end();
        second.generate(start, i, 2, b, P, ss);
        newL.push_back(second);
    }
    newL.sort();
    vector<circuit> V3(newL.begin(), newL.end());
    clock_t e3 = clock();
    cout << (e3 - s3) << endl;

    // fourth situation
    vector<circuit> V4;
    clock_t s4 = clock();
    for (UINT i = 0; i < 0xfff; ++i)
    {
        circuit second;
        bool b = find(Y.begin(), Y.end(), i) != Y.end();
        second.generate(start, i, 2, b, P, ss);
        V4.push_back(second);
    }
    sort(V4.begin(), V4.end());
    clock_t e4 = clock();
    cout << (e4 - s4) << endl;

    for (auto ele : V2)
        result_file  << ele;
    */
    result_file.close();

    return 0;
}

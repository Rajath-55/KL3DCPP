from typing import Dict, List, Optional, Union

from header import *

# Header definitions
TILE_WIDTH = 0.0025  # IN DIRECTION OF COLUMN
TILE_HEIGHT = 0.0025  # IN DIRECTION OF ROWS
NO_OF_RUN = 5
NUM_LAYERS = 2
# no_particle = 100


class GraphInfo:
    def __init__(self) -> None:
        self.No_nodes = 0
        self.actual_no_of_nodes = 0
        self.rows = 0
        self.colums = 0
        self.edges = 0


graph_info = GraphInfo()

# check this part
cost_best = 99999999
a = [[]]
save_partition = [[]]
seed_value = []

"""
Function Read_graph() 
Inputs: name => char array containing the name of the input graph_info 
        dist => control variable telling if the placement should be Distributive or Non-Distributive 
Output: netlist => 2D adjacency matrix holding float values 
        
Description : Function reads the input graph and saves the bandwidth requirements of cores as a 2D adjacency matrix (netlist). 
Since K-L is a bipartitioning algorithm, in the graphs where the no. of nodes is not a power of 2, dummy nodes are added. The communication bandwidth between these dummy nodes and the actual nodes is set to be zero, while the communication bandwidth between dummy nodes depends on the user. The user is given 2 choices in the type of placement he wants namely, Distributed(1) and Non-Distributed(0). For Distributed type placement the communication bandwidth between the dummy nodes is taken to be zero, while for Non-Distributed type placement the communication bandwith between dummy nodes is taken to be INFINITE, ie. a very large number defined as "MAX". The function also saves the different attributes of the graph in a global structure. 
"""


def Read_graph(name: str, dist: str) -> Union[List[List[float]], None]:
    bandwidth_graph_pointer = open(name, "r+")
    nodes = int(bandwidth_graph_pointer.readline().strip())
    graph_info.actual_no_of_nodes = nodes
    actual_no_of_nodes = nodes

    check_power_of_2_nodes = math.log(actual_no_of_nodes) / math.log(2.0)
    power_of_2_nearest = int(check_power_of_2_nodes)

    logging.info(
        f"Check number of nodes {check_power_of_2_nodes} Nearest power of 2 {power_of_2_nearest}"
    )

    if power_of_2_nearest > float(check_power_of_2_nodes):
        logging.info("Adding dummy nodes.")
        new_power_of_two = pow(2, check_power_of_2_nodes + 1)
        new_nodes = int(new_power_of_two)
        logging.info(f"New Nodes with dummy - {new_nodes}")

    netlist = [[0.0 for _ in range(new_nodes)] for _ in range(new_nodes)]

    for row in range(actual_no_of_nodes):
        buff = bandwidth_graph_pointer.readline().strip().split()

        for column, value in enumerate(buff):
            if value and "INF" not in value:
                z = float(value)
            else:
                z = 0.0
            netlist[row][column] = z

    temp_pow = math.log(nodes) / math.log(2.0)
    split1 = math.floor(temp_pow / 2.0)
    split2 = math.ceil(temp_pow / 2.0)
    graph_info.No_nodes = nodes
    graph_info.rows = int(pow(2, split1))
    graph_info.colums = int(pow(2, split2))

    logging.info(f"rows {graph_info.rows}, cols {graph_info.colums}")
    bandwidth_graph_pointer.close()
    distri = int(dist)

    if power_of_2_nearest > float(check_power_of_2_nodes) and distri == 1:
        # Fill 0s in dummy rows
        for i in range(actual_no_of_nodes, nodes):
            for j in range(nodes):
                if j >= actual_no_of_nodes:
                    netlist[i][j] = 0.0
                else:
                    netlist[i][j] = 0.0

        # Fill 0s in dummy cols
        for i in range(actual_no_of_nodes):
            for j in range(actual_no_of_nodes, nodes):
                netlist[i][j] = 0.0

        # Fill 0s in Diagonal
        for i in range(nodes):
            netlist[i][i] = 0.0

    elif power_of_2_nearest > float(check_power_of_2_nodes):
        for i in range(actual_no_of_nodes, nodes):
            for j in range(nodes):
                if j >= actual_no_of_nodes:
                    netlist[i][j] = MAX
                else:
                    netlist[i][j] = 0.0

        for i in range(actual_no_of_nodes):
            for j in range(actual_no_of_nodes, nodes):
                netlist[i][j] = 0.0

        for i in range(nodes):
            netlist[i][i] = 0.0

    return netlist


"""
Function initialize_add() 
Inputs: NO INPUTS 
Output: NO OUTPUTS 
        
Description : This initializes the Global addresses for a Mesh topology for given number of nodes. It also generates the partition id reference table according to the K-L algorithm. Thus, the K-L algorithm uses this partition ID reference table and does not generate the partition ids at each iteration saving computaional time.  
"""


def initialize_add() -> None:
    nodes = graph_info.No_nodes
    global a
    ht_add = []
    row_add = []
    col_add = []
    log_of_nodes = math.log(float(nodes)) / math.log(2.0)
    n = int(log_of_nodes)
    temp_add = [[0 for _ in range(n - 1)] for _ in range(n - 1)]

    for i in range(1, n):
        c = 0
        l = 0
        for j in range(2**i):
            for k in range(2 ** (n - i)):
                temp_add[i - 1][l] = c
                l += 1
            c += 1
    for i in range(0, nodes, 4):
        row_add[i] = col_add[i] = row_add[i + 1] = col_add[i + 2] = 0
        col_add[i + 1] = row_add[i + 2] = col_add[i + 3] = row_add[i + 3] = 1

    count = 0
    xyz = 1
    for i in range(4, nodes, i * 4):
        count = 0
        xyz = xyz * 2
        j = 0
        while i * count < nodes:
            j = j % 4
            for k in range(i):
                if k < i and i * count + k < nodes:
                    if j == 2 or j == 3:
                        row_add[i * count + k] += xyz
                    if j == 1 or j == 3:
                        col_add[i * count + k] += xyz

            j += 1
            count += 1

    c = 0
    logging.info(f"{nodes//NUM_LAYERS}")
    for i in range(nodes):
        if i % (nodes // NUM_LAYERS) == 0:
            c += 1
        ht_add[i] = c - 1

    address.colum = col_add
    address.row = row_add
    address.height = ht_add

    logging.info("Row addressing - " + str(address.row))
    logging.info("Column addressing - " + str(address.column))
    logging.info("Height addressing - " + str(address.height))
    a = temp_add


"""
Function main() 
Inputs: argc => number of command line inputs 
          argv => array of strings ie. 2D array of charactres, it contains the command line inputs. argv[1] holds the name of input graph file, argv[2] holds the initial number of cuts, argv[3] holds if the placement is to be Distributive(1) or Non-Distributive(0) 
Output: NO OUTPUTS 
        
Description : Topmost function. Controls all the other functions and saves the results
 
NOTE : FOR RUNNING THE PROGRAM THE INPUT GRAPHS MUST BE IN THE SAME FOLDER AS THE MAIN AND HEADER FILES. THE WAY OF COMMAND LINE SEQUENCE FOR EXECUTION IS AS FOLLOWS:
./a.out <Graph name> <initial number of cuts> <Distributive(1)/Non-Distributive(0)>
ex.: 	./a.out Graph4.txt 100 0  		
"""


def main(argc: int, argv: List[str]) -> int:
    global NUM_LAYERS
    curr_time = time.process_time()

    if len(argv) < 2:
        logging.critical(
            "Use: python3 main_comd1.py input_graph.txt no_of_runs is_distributive num_layers"
        )
        exit(0)
    pass


def perform_KL(
    argv: List[str], netlist: List[List[float]], output: List[float]
) -> None:
    """
    Function perform_KL() 
    Inputs: argv[] => holds the command line inputs 
            netlist => adjacency matrix 
            output => Since there are multiple outputs this array is passed by reference from the controlling function to save results 
            Core_Power => float array holding the power values of the cores 
    Output: output => float array holding results generated by the algorithm. Passed by reference from controlling function 
            
    Description : This function is the controlling function of all the partitioning, mapping, cost and thermal (hotspot) functions. Since the results from K-L algorithm depend on initial cut, a number of initial random cuts are generated and then K-L algorithm is called to generate the result for each initial cut. Each of these result is saved and the best result is found out. The best result is returned to the main() function.
    """
    pass


def map_nodes(
    nodes: int, final_partition_core: List[int], graph: List[List[float]]
) -> float:
    """
    Function map_nodes()
    Inputs: nodes => number of nodes in the graph
            final_partition_core => array holding the final core sequence generated after partitioning
            graph => adjacency matrix
            core_pow => array containing power values of each core
    Output: final communication cost of the mapping generated
            final_partition_core => It holds the final core sequence generated after different cost and power improvements

    Description : This function takes the core sequence generated from the partioning algorithm and first applies iterative power improvements (local followed by global) and the iterative cost improvements (local followed by global). Then the final best cost thus found and the corresponding core sequence is returned.
    """
    temp_final_partition_core = [0 for _ in range(nodes)]
    best_cost = [0.0 for _ in range(4)]
    ############## Iterative improvement phase for Communication cost improvements ################

    # Local phase

    iterative_improvement(
        graph, final_partition_core, nodes, 1
    )  # Local swap analysis. '1' stands for local.

    # global phase

    Final_Global_best_cost = cost(final_partition_core, graph, nodes)

    for i in nodes:
        temp_final_partition_core[i] = final_partition_core[i]

    iterative_improvement(graph, temp_final_partition_core, nodes, 0)
    best_cost[3] = cost(temp_final_partition_core, graph, nodes)
    best_cost[1] = best_cost[2] = best_cost[0] = MAX

    while 1:
        iterative_improvement(graph, temp_final_partition_core, nodes, 0)
        temp = cost(temp_final_partition_core, graph, nodes)
        print(f"\n Cost : {best_cost[3]}")

        if temp <= best_cost[3]:
            best_cost[0] = best_cost[1]
            best_cost[1] = best_cost[2]
            best_cost[2] = best_cost[3]
            best_cost[3] = temp
            # The Global analysis continues till there is no improvement in cost for 4 consecutive passes.
            # The number of comparisions for communication cost can be increased if required.

        if (
            temp == best_cost[3]
            and temp == best_cost[2]
            and temp == best_cost[1]
            and temp == best_cost[0]
        ):
            break

    if Final_Global_best_cost > best_cost[3]:
        Final_Global_best_cost = best_cost[3]
        for i in nodes:
            final_partition_core[i] = temp_final_partition_core[i]

    Final_Global_best_cost = cost(final_partition_core, graph, nodes)

    del temp_final_partition_core
    return Final_Global_best_cost



def iterative_improvement(
    graph: List[List[float]], final_partition_core: List[int], nodes: int, local: int
) -> None:
    """
    Function iterative_improvement() 
    Inputs: nodes => number of nodes in the graph 
            final_partition_core => array holding the final core sequence generated after partitioning 
            graph => adjacency matrix 
            local => control varaible which controls the type of pass local or global 
    Output: final_partition_core => It holds the final core sequence generated after cost improvement 
            
    Description : This function takes the core sequence generated from the partitioning and performs iterations for improving communication cost of the mapping. The partitions are selected in pairs, at a level of partitioning, and different arrangements are generated in one of the partition keeping the other fixed. The same procedure is repeated for all the levels of partitioning. 
    """

    n4 = math.log2(nodes)
    n = int(n4)
    level = int(math.log2(nodes)) -1

    temp_final_partition_core = [[0 for _ in range(nodes)] for _ in range(4)]
    Global_best = 0.0
    best_cost = [0.0 for _ in range(4)]
    temp_row = [0 for _ in range(nodes)]
    temp_col = [0 for _ in range(nodes)]
    temp_ht = [0 for _ in range(nodes)]

    for curr_lvl in range(level,0,-1):
        for i in range(nodes):
            temp_final_partition_core[0][i] = final_partition_core[i]
            temp_final_partition_core[1][i] = final_partition_core[i]
            temp_final_partition_core[2][i] = final_partition_core[i]
            temp_final_partition_core[3][i] = final_partition_core[i]

        for i in range(int(math.pow(2, curr_lvl)),2):
            for j in range(nodes):
                temp_row[j] = address.row[j]
                temp_col[j] = address.column[j]
                temp_ht[j] = address.height[j]
            
            for k in range(nodes):
                if a[curr_lvl-1][k] == i+1:
                    break
            
            t = int(math.pow(2, n - curr_lvl))

            min_row = address.row[k-t]
            min_col = address.column[k-t]
            max_row = address.row[k-t]
            max_col = address.column[k-t]
            max_ht = address.height[k-t]
            min_ht = address.height[k-t]

            for w in range(k-t, k):
                if address.row[w] > max_row:
                    max_row = address.row[w]
                if address.column[w] > max_col:
                    max_col = address.column[w]
                if address.height[w] > max_ht:
                    max_ht = address.height[w]

            avg_row = (max_row + min_row) / 2.0
            avg_col = (max_col + min_col) / 2.0
            avg_ht = (max_ht + min_ht) / 2.0

            best_cost[0] = flip(graph, temp_final_partition_core[0], k, t, nodes, local)

            if curr_lvl == level:
                #### Flip along lvl axis ####
                difference = (max_ht - min_ht)/2.0
                shift = math.ceil(difference)

                for w in range(k-t, k):
                    if address.height[w] >= avg_ht:
                        address.height[w] -= shift
                    elif address.height[w] < avg_ht:
                        address.height[w] += shift
                
                for w in range(nodes):
                    for j in range(nodes):
                        if temp_row[w]==address.row[j] and temp_col[w]==address.colum[j] and temp_ht[w]==address.height[j] :
                            temp_final_partition_core[1][w] = final_partition_core[j]
                            break
                
                for j in range(nodes):
                    address.row[j] = temp_row[j]
                    address.column[j] = temp_col[j]
                    address.height[j] = temp_ht[j]

                best_cost[1] = flipd(graph, temp_final_partition_core[1], k, t, nodes, local)


                ##### Flip along horizontal axis ####

                difference = (max_row - min_row)/2.0
                shift = math.ceil(difference)

                for w in range(k-t, k):
                    if address.row[w] >= avg_row:
                        address.row[w] -= shift
                    elif address.row[w] < avg_row:
                        address.row[w] += shift

                for w in range(nodes):
                    for j in range(nodes):
                        if temp_row[w]==address.row[j] and temp_col[w]==address.colum[j] and temp_ht[w]==address.height[j] :
                            temp_final_partition_core[2][w] = final_partition_core[j]
                            break

                for j in range(nodes):
                    address.row[j] = temp_row[j]
                    address.column[j] = temp_col[j]
                    address.height[j] = temp_ht[j]

                best_cost[2] = flipd(graph, temp_final_partition_core[2], k, t, nodes, local)


                ##### Flip along lvl axis ####

                difference = (max_ht - min_ht)/2.0
                shift = math.ceil(difference)

                for w in range(k-t, k):
                    if address.height[w] >= avg_ht:
                        address.height[w] -= shift
                    elif address.height[w] < avg_ht:
                        address.height[w] += shift

                for w in range(nodes):
                    for j in range(nodes):
                        if temp_row[w]==address.row[j] and temp_col[w]==address.colum[j] and temp_ht[w]==address.height[j] :
                            temp_final_partition_core[3][w] = final_partition_core[j]
                            break

                for j in range(nodes):
                    address.row[j] = temp_row[j]
                    address.column[j] = temp_col[j]
                    address.height[j] = temp_ht[j]

                best_cost[3] = flipd(graph, temp_final_partition_core[3], k, t, nodes, local)
            
            else:
                #### Flip along horizontal axis ####

                difference = (max_row - min_row)/2.0
                shift = math.ceil(difference)

                for w in range(k-t, k):
                    if address.row[w] >= avg_row:
                        address.row[w] -= shift
                    elif address.row[w] < avg_row:
                        address.row[w] += shift

                for w in range(nodes):
                    for j in range(nodes):
                        if temp_row[w]==address.row[j] and temp_col[w]==address.colum[j] and temp_ht[w]==address.height[j] :
                            temp_final_partition_core[1][w] = final_partition_core[j]
                            break


                for j in range(nodes):
                    address.row[j] = temp_row[j]
                    address.column[j] = temp_col[j]
                    address.height[j] = temp_ht[j]

                best_cost[1] = flipd(graph, temp_final_partition_core[1], k, t, nodes, local)

                ##### Flip along vertical axis #### 

                difference = (max_col - min_col)/2.0
                shift = math.ceil(difference)

                for w in range(k-t, k):
                    if address.column[w] >= avg_col:
                        address.column[w] -= shift
                    elif address.column[w] < avg_col:
                        address.column[w] += shift

                for w in range(nodes):
                    for j in range(nodes):
                        if temp_row[w]==address.row[j] and temp_col[w]==address.colum[j] and temp_ht[w]==address.height[j] :
                            temp_final_partition_core[2][w] = final_partition_core[j]
                            break
                
                for j in range(nodes):
                    address.row[j] = temp_row[j]
                    address.column[j] = temp_col[j]
                    address.height[j] = temp_ht[j]

                best_cost[2] = flipd(graph, temp_final_partition_core[2], k, t, nodes, local)

                ##### Flip along horizontal axis ####

                difference = (max_row - min_row)/2.0
                shift = math.ceil(difference)

                for w in range(k-t, k):
                    if address.row[w] >= avg_row:
                        address.row[w] -= shift
                    elif address.row[w] < avg_row:
                        address.row[w] += shift

                for w in range(nodes):
                    for j in range(nodes):
                        if temp_row[w]==address.row[j] and temp_col[w]==address.colum[j] and temp_ht[w]==address.height[j] :
                            temp_final_partition_core[3][w] = final_partition_core[j]
                            break
                
                for j in range(nodes):
                    address.row[j] = temp_row[j]
                    address.column[j] = temp_col[j]
                    address.height[j] = temp_ht[j]

                best_cost[3] = flipd(graph, temp_final_partition_core[3], k, t, nodes, local)

            ###### Checked all combinations ########
                
                Global_best = best_cost[0]
                best_partition = 0
                for j in range(4):
                    if best_cost[j] < Global_best:
                        Global_best = best_cost[j]
                        best_partition = j
                
                for j in range(nodes):
                    final_partition_core[j] = temp_final_partition_core[best_partition][j]

    del temp_col
    del temp_row
    del temp_final_partition_core        


def flip(
    G: List[List[float]],
    final_partition_core: List[int],
    k: int,
    t: int,
    nodes: int,
    local: int,
) -> float:
    """
    Function flip() 
    Inputs: nodes => number of nodes in the graph 
            final_partition_core => array holding the final core sequence generated after partitioning 
            G => adjacency matrix 
            local => control varaible which controls the type of pass local or global 
            k => start index of the partition under consideration
            t => size of the partition under consideration
    Output: final_partition_core => It holds the final core sequence generated after cost improvement 
            
    Description : This function takes a partition and rearranges the core within it and checks for cost improvements. The arrangement with the lowest cost is returned.
    """
    return 0.0


def flipd(
    G: List[List[float]],
    final_partition_core: List[int],
    k: int,
    t: int,
    nodes: int,
    local: int,
) -> float:
    return 0.0

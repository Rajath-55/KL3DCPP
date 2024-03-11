import random
from typing import List, Dict, Union, Optional
import math
import time
import logging, sys 
import numpy as np

# Logging info
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(funcName)s %(lineno)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    filename="info.log",
    filemode = 'w',
    encoding="utf-8",
    level=logging.DEBUG,
)

# CAVEAT: q directly is not being updated. Using a class for this purpose.
class new_q:

    def __init__(self):
        self.q = 0
    
    def update(self, val):
        self.q = val
    
    def get(self):
        return self.q
    

q = new_q()

# Header definitions
MAX = 1000000
LAMBDA = 2
WEIGHT = 1.0
header_q = 0
seed_value = [0]*1000
save_partition = None


class VirtualAddress:
    def __init__(self) -> None:
        self.row = []
        self.column = []
        self.height = []


address = VirtualAddress()

# floating point initializations
init_comm_cut_cost = 0
init_powbal = 0
wt = WEIGHT
wt1 = WEIGHT


def KL_partition(
    index: int,
    graph: List[List[float]],
    n: int,
    core_id: List[int],
    final_partition_core: List[int],
) -> None:
    """
    Function KL_partition()
    Inputs: n => number of nodes in the graph
            final_partition_core => array holding the final core sequence generated after partitioning
            graph => adjacency matrix
            core_id => array containing the initial sequence of cores
            core_pow => array containing power values of each core

    Output: final_partition_core => array holding the final core sequence generated after partitioning
    Description : This function is the control function for KL() partitioning function. It recursively calls itself till the smallest partition size is reached. The smallest size of partition is controlled by 'lambda' which is '2' for MoT and Mesh and '4' for BFT. It takes a number of random cuts for each partition by calling KL() and calculates the cut cost. The partitioning with the best cut cost is accepted and passed for further partitioning
    """
    best_cost, cut = sys.maxsize, 0.0
    global header_q

    if n <= LAMBDA:
        # Lambda = 2 for MoT and Mesh, while Lambda = 4 for BFT
        # Lambda represents the smallest partitions size.
        logging.info("LAMBDA = " + str(LAMBDA) + " core_id size = " + str(len(core_id)))
        for i in range(LAMBDA):
            # print("q = " +str(q.get()))
            final_partition_core[index][q.get()] = core_id[i]
            q.update(q.get() + 1)
        return

    partition = [[0 for _ in range(n // 2)], [0 for _ in range(n // 2)]]
    best_partition = [[0 for _ in range(n // 2)], [0 for _ in range(n // 2)]]

    for i in range(
        n // 4
    ):  # number of random cuts reduced to n/4. Can be changed according to requirement
        KL(index, graph, core_id, n, partition)

        cut = partition_cost(graph, partition[0], partition[1], n // 2)

        if cut < best_cost:
            best_cost = cut
            for j in range(n // 2):
                best_partition[0][j] = partition[0][j]
                best_partition[1][j] = partition[1][j]

    del partition
    logging.info(f"best cost={best_cost}")

    KL_partition(index, graph, n // 2, best_partition[0], final_partition_core)
    KL_partition(index, graph, n // 2, best_partition[1], final_partition_core)


def random_partition(
    core_id: List[int], n: int, indx: int, partition: List[List[int]]
) -> None:
    """
    Function random_partition()
    Inputs: n => number of nodes in the graph
            partition => 2D array holding the 2 partitions generated
            core_id => array containing the initial sequence of cores

    Output: partition => 2D array holding the 2 partitions generated

    Description : This function randomly bi-partitions the core sequence passed to it.
    """
    j = k = 0
    random.seed(seed_value[indx])
    
    for i in range(n):
        rand = random.randint(0, 5147483647)
        
        if rand%2 == 0:
            if j < n//2:
                partition[0][j] = core_id[i]
                j += 1
            else:
                partition[1][k] = core_id[i]
                k += 1
        else:
            if k < n//2:
                partition[1][k] = core_id[i]
                k += 1
            else:
                partition[0][j] = core_id[i]
                j += 1
                
        


def KL(
    indx: int,
    graph: List[List[float]],
    core_id: List[int],
    n: int,
    partition: List[List[int]],
) -> None:
    """
    Function KL()
    Inputs: n => number of nodes in the graph
            partition => 2D array holding the 2 partitions generated by KL algorithm
            graph => adjacency matrix
            core_id => array containing the initial sequence of cores
            core_pow => array containing power values of each core

    Output: partition => 2D array holding the 2 partitions generated by KL algorithm

    Description : This function is the implementation of the KL bi-patitioning algorithm. Firstly, it calls random_partition() to generate 2 random partitions and then performs swapping according to KL algorithm. The partitions thus generated are returned.
    """
    random_partition(core_id, n, indx, partition)
    
    
    d_a = d_b = [0.0] * (n//2)
    gain_k = sum_gain_k = [0.0] * (n//2)
    
    temp_a = [partition[0][i] for i in range(n//2)]
    temp_b = [partition[1][i] for i in range(n//2)]
    
    counter = flag = 0
    old_maximum = float('inf')
    
    global init_comm_cut_cost
    
    if init_comm_cut_cost == 0:
        init_comm_cut_cost = partition_cost(graph, temp_a, temp_b, n // 2)
        
    while True:
        counter = 0
        
        while counter < n//2 : 
            for j in range(counter, n//2):
                int_cost_a = ext_cost_a = int_cost_b = ext_cost_b = 0.0
                
                for i in range(n//2):
                    int_cost_a += graph[temp_a[j]][temp_a[i]]
                    ext_cost_a += graph[temp_a[j]][temp_b[i]]
                    
                    int_cost_b += graph[temp_b[j]][temp_b[i]]
                    ext_cost_b += graph[temp_b[j]][temp_a[i]]
                    
                d_a[j] = ext_cost_a - int_cost_a
                d_b[j] = ext_cost_b - int_cost_b

            gain = -float('inf')
                
            for i in range(counter, n//2):
                for j in range(counter, n//2):
                    temp = (d_a[i] + d_b[j] - graph[temp_a[i]][temp_b[j]] - graph[temp_b[j]][temp_a[i]])/init_comm_cut_cost
                    if temp > gain:
                        gain = temp
                        new_i, new_j = i, j
                    
                gain_k[counter] = gain
                counter += 1
                swap_a, swap_b = temp_a[new_i], temp_b[new_j]
                    
                for i in range(new_i, counter - 1, -1):
                    temp_a[i] = temp_a[i - 1]
                
                temp_a[counter - 1] = swap_b
                    
                for i in range(new_j, counter - 1, -1):
                    temp_b[i] = temp_b[i - 1]
                    
                temp_b[counter - 1] = swap_a
                
        for i in range(counter):
            sum_gain_k[i] = 0.0
            for j in range(i, -1, -1):
                sum_gain_k[i] += gain_k[j]
                
        maxi = sum_gain_k[0]
        k = 0
        for i in range(counter):
            if sum_gain_k[i] > maxi:
                maxi = sum_gain_k[i]
                k = i

        for i in range(k + 1):
            for j in range(n // 2):
                if temp_b[i] == partition[0][j]:
                    partition[0][j] = temp_a[i]

                if temp_a[i] == partition[1][j]:
                    partition[1][j] = temp_b[i]
                
        temp_a =[partition[0][i] for i in range(n//2)]
        temp_b = [partition[1][i] for i in range(n//2)]
                
        if flag == 1:
            return
                
        if old_maximum <= maxi:
            flag = 1
                
        old_maximum = maxi
                
        if maxi < 0 or k >= (n//2 - 1):
            break
        
    return
                        
                        
            
def partition_cost(
    graph: List[List[float]], A: List[int], B: List[int], n: int
) -> float:
    """
    Function partition_cost()
    Inputs: n => number of nodes in the graph
            graph => adjacency matrix
            A => holds the cores present in partition A
            B => holds the cores present in partition B

    Output: cut cost or the cost of inter-communication between partition A and B

    Description : This function calculates the cut cost between partition A and B
    """
    cut = 0.0

    for i in range(n):
        for j in range(n):
            cut += graph[A[i]][B[j]]

    # print("Cut= ", cut)
    return cut


# FOR MESH (Mesh router related code)


def cost_local(
    map: List[int], G: List[List[float]], n: int, start: int, end: int
) -> float:
    """
    Function cost_local()
    Inputs: n => number of nodes in the graph
            G => adjacency matrix
            start => starting index of partition under consideration
            end => end index of the partition under consideration
            map => suquence of cores
    Output: local communication cost (in Hops x Bandwidth)

    Description : This function calculates the local communication cost (in Hops x Bandwidth) of a partition
    """
    cost = 0.0

    for index1 in range(start, end):
        for index2 in range(index1 + 1, end):
            if G[map[index1]][map[index2]] != MAX and G[map[index1]][map[index2]] != 0:
                if G[map[index1]][map[index2]] == G[map[index2]][map[index1]]:
                    cost += float(hops(index1, index2, n) * G[map[index1]][map[index2]])
                else:
                    cost += float(
                        hops(index1, index2, n) * G[map[index1]][map[index2]]
                    ) + float(hops(index1, index2, n) * G[map[index2]][map[index1]])

    return cost


def cost(map: List[int], G: List[List[float]], n: int) -> float:
    """
    Function cost()
    Inputs: n => number of nodes in the graph
            G => adjacency matrix
            map => suquence of cores
    Output: communication cost (in Hops x Bandwidth) of the whole network

    Description : This function calculates the total communication cost (in Hops x Bandwidth) of the network
    """
    cost = 0.0

    for index1 in range(n):
        for index2 in range(index1 + 1, n):
            if G[map[index1]][map[index2]] != MAX and G[map[index1]][map[index2]] != 0:
                if G[map[index1]][map[index2]] == G[map[index2]][map[index1]]:
                    cost += float(hops(index1, index2, n) * G[map[index1]][map[index2]])
                else:
                    cost += float(
                        hops(index1, index2, n) * G[map[index1]][map[index2]]
                    ) + float(hops(index2, index1, n) * G[map[index2]][map[index1]])

    return cost


def cost_cc(map: List[int], G: List[List[float]], n: int) -> float:
    """
    Function cost_cc()
    Inputs: n => number of nodes in the graph
            G => adjacency matrix
            map => suquence of cores
    Output: communication cost (in Router cycles x Bandwidth) of the whole network

    Description : This function calculates the total communication cost (in Router cycles x Bandwidth) of the network
    """
    cost = 0.0

    for index1 in range(n):
        for index2 in range(index1 + 1, n):
            if G[map[index1]][map[index2]] != MAX and G[map[index1]][map[index2]] != 0:
                if G[map[index1]][map[index2]] == G[map[index2]][map[index1]]:
                    cost += float(
                        (hops(index1, index2, n) + 1) * 2 * G[map[index1]][map[index2]]
                    )
                else:
                    # Router cycles are taken as 2
                    cost += float(
                        (hops(index1, index2, n) + 1) * 2 * G[map[index1]][map[index2]]
                    ) + float(
                        (hops(index1, index2, n) + 1) * 2 * G[map[index2]][map[index1]]
                    )

    return cost


def hops(snode: int, dnode: int, n: int) -> int:
    """
    Function Hops()
    Inputs: n => number of nodes in the graph
            node1 => starting node
            node2 => Destination node
    Output: hop_count => noumber of hops between node1 and node2

    Description : This function calculates the total Hops required to reach node2 starting from node1
    """
    hop_count = 0
    srow = address.row[snode]
    scol = address.column[snode]
    drow = address.row[dnode]
    dcol = address.column[dnode]
    sheight = address.height[snode]
    dheight = address.height[dnode]

    hop_count = abs(sheight - dheight) + abs(srow - drow) + abs(scol - dcol)
    return hop_count


def diff(first: int, second: int) -> float:
    """
    Function diff()
    Inputs: first => integer
            second => integer
    Output: floating point absolute difference of the two integers.

    Description: Wrapper functionn to get diff between two integers
    """
    return float(abs(first - second))

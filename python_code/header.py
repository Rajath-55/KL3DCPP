from typing import List, Dict, Union, Optional

# Header definitions
MAX = 1000000
LAMBDA = 2
WEIGHT = 1.0
header_q = 0

class VirtualAddress:
    def __init__(self) -> None:
        self.row = []
        self.column = []
        self.height = []

address = VirtualAddress()

#floating point initializations
init_comm_cut_cost = 0
init_powbal = 0
wt = WEIGHT
wt1 = WEIGHT

'''
Function KL_partition()
Inputs: n => number of nodes in the graph 
		final_partition_core => array holding the final core sequence generated after partitioning 
		graph => adjacency matrix 
		core_id => array containing the initial sequence of cores 
		core_pow => array containing power values of each core 

Output: final_partition_core => array holding the final core sequence generated after partitioning 
Description : This function is the control function for KL() partitioning function. It recursively calls itself till the smallest partition size is reached. The smallest size of partition is controlled by 'lambda' which is '2' for MoT and Mesh and '4' for BFT. It takes a number of random cuts for each partition by calling KL() and calculates the cut cost. The partitioning with the best cut cost is accepted and passed for further partitioning
'''
def KL_partition(index : int, graph : List[List[float]], n : int, core_id : List[int], final_partition_core : List[int]) -> None:
    pass

'''
Function random_partition() 
Inputs: n => number of nodes in the graph 
		  partition => 2D array holding the 2 partitions generated 
		  core_id => array containing the initial sequence of cores 

Output: partition => 2D array holding the 2 partitions generated 
		
Description : This function randomly bi-partitions the core sequence passed to it.
'''
def random_partition(core_id : List[int], n : int, indx : int, partition : List[List[int]]) -> None:
    pass

'''
 Function KL()
Inputs: n => number of nodes in the graph
		  partition => 2D array holding the 2 partitions generated by KL algorithm 
		  graph => adjacency matrix 
		  core_id => array containing the initial sequence of cores
		  core_pow => array containing power values of each core 

Output: partition => 2D array holding the 2 partitions generated by KL algorithm 
		
Description : This function is the implementation of the KL bi-patitioning algorithm. Firstly, it calls random_partition() to generate 2 random partitions and then performs swapping according to KL algorithm. The partitions thus generated are returned.
'''
def KL(indx : int, graph : List[List[float]], core_id : List[int], n : int, partition : List[List[int]]) -> None:
    pass


'''
Function partition_cost() 
Inputs: n => number of nodes in the graph 
		  graph => adjacency matrix 
		  A => holds the cores present in partition A 
		  B => holds the cores present in partition B 

Output: cut cost or the cost of inter-communication between partition A and B 
		
 Description : This function calculates the cut cost between partition A and B
'''
def partition_cost(graph : List[List[float]], A : List[int], B : List[int], n : int) -> float:
    pass


# FOR MESH (Mesh router related code)


'''
 Function cost_local() 
Inputs: n => number of nodes in the graph 
		  G => adjacency matrix 
		  start => starting index of partition under consideration
		  end => end index of the partition under consideration 
		  map => suquence of cores
Output: local communication cost (in Hops x Bandwidth) 
		
 Description : This function calculates the local communication cost (in Hops x Bandwidth) of a partition
'''
def cost_local(map : List[int], G : List[List[float]], n : int, start : int, end : int) -> float:
    pass


'''
 Function cost() 
Inputs: n => number of nodes in the graph 
		  G => adjacency matrix 
		  map => suquence of cores
Output: communication cost (in Hops x Bandwidth) of the whole network
		
 Description : This function calculates the total communication cost (in Hops x Bandwidth) of the network 
'''
def cost(map : List[int], G : List[List[float]], n : int) -> float:
    pass

'''
 Function cost_cc() 
Inputs: n => number of nodes in the graph 
		  G => adjacency matrix 
		  map => suquence of cores
Output: communication cost (in Router cycles x Bandwidth) of the whole network
		
 Description : This function calculates the total communication cost (in Router cycles x Bandwidth) of the network 
'''
def cost_cc(map : List[int], G : List[List[float]], n : int) -> float:
    pass


'''
 Function Hops() 
Inputs: n => number of nodes in the graph 
		  node1 => starting node
		  node2 => Destination node
Output: hop_count => noumber of hops between node1 and node2
		
 Description : This function calculates the total Hops required to reach node2 starting from node1 
'''
def Hops(snode : int, dnode : int, n : int) -> int:
    hop_count = 0
    srow = address.row[snode]
    scol = address.column[snode]
    drow = address.row[dnode]
    dcol = address.column[dnode]
    
    if address.height[snode] != address.ht[dnode]:
        hop_count += 1
    hop_count = hop_count + abs(srow - drow) + abs(scol - dcol)
    return hop_count
    
'''
Function diff()
Inputs: first => integer 
		second => integer
Output: floating point absolute difference of the two integers.

Description: Wrapper functionn to get diff between two integers
'''
def diff(first : int, second : int) -> float:
    return float(abs(first-second))
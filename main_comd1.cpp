 #include "Header.h"
 #include <iomanip>

#define TILE_WIDTH 0.0025	//IN DIRECTION OF COLUMN
#define TILE_HEIGHT 0.0025	//IN DIRECTION OF ROWS
#define NO_OF_RUN 5
//#define no_particle 100

/**********Global variables***********/
int **a;
float cost_best=99999999;
int **save_partition;
int *seed_value;
int NUM_LAYERS = 1;
struct graph_info
{
    unsigned short No_nodes;
    unsigned short actual_No_nodes;
    unsigned short rows;
    unsigned short colums;
    unsigned short edges;

}Graph_inf;




/************************************************ Declaration of the functions ************************************************/
//void iteration_power(float **graph, int *final_partition_core,int nodes, int local, float *core_pow);
//float variance(float *core_pow, int *final_partition_core, int start, int end);
void perform_KL(char *argv[], float **netlist, float * output);
float map_nodes(int nodes, int *final_partition_core, float **graph);
float flipd(float **G,int *final_partition_core,int k,int t,int nodes,int local);
/*****************************************************************************************************************************/

/*Function to print array*/
void print_arr(int *arr, int len){
	printf("[ ");
	for(int i = 0; i < len; ++i){
		printf("%d ", arr[i]);
	}
	printf("]\n");
}

/* Function Read_graph() */
/*Inputs: name => char array containing the name of the input graph_info */
/*        dist => control variable telling if the placement should be Distributive or Non-Distributive */
/*Output: netlist => 2D adjacency matrix holding float values */
		
/* Description : Function reads the input graph and saves the bandwidth requirements of cores as a 2D adjacency matrix (netlist). 
 Since K-L is a bipartitioning algorithm, in the graphs where the no. of nodes is not a power of 2, dummy nodes are added. The communication bandwidth between these dummy nodes and the actual nodes is set to be zero, while the communication bandwidth between dummy nodes depends on the user. The user is given 2 choices in the type of placement he wants namely, Distributed(1) and Non-Distributed(0). For Distributed type placement the communication bandwidth between the dummy nodes is taken to be zero, while for Non-Distributed type placement the communication bandwith between dummy nodes is taken to be INFINITE, ie. a very large number defined as "MAX". The function also saves the different attributes of the graph in a global structure. */
 
float ** Read_graph(char *name, char *dist)
{
	float **netlist;
	int i,j,row,co,count,n1,n4,nodes,actual_no_of_nodes;
	float z,n2,n3;
	char buff[5];
	fstream fp1;
	fp1.open(name,ios::in);
	

    fp1>>nodes;
    fp1>>buff;
   
	Graph_inf.actual_No_nodes= actual_no_of_nodes=nodes;
	Graph_inf.edges=0;

	n2=log((float)nodes)/log(2.0);
	n1=(int)n2;
	printf("%d %f %d",n1,n2,nodes);
		if(n2>(float)n1)
		{
			printf("\n\t\tAdding dummy nodes\n");
			n3=pow(2.0,n1+1);
			nodes=(int)n3;
			printf("\nn=%d",nodes);
		}

	row=0;
	co=0;
	count=0;
	netlist=(float **)malloc(nodes*sizeof(float *));
		for(i=0;i<nodes;i++)
			netlist[i]=(float *)malloc(nodes*sizeof(float));

		while(count<actual_no_of_nodes*actual_no_of_nodes)
		{
		
	   		z=atof(buff);
			if(!strcmp(buff,"INF"))
				z=0;
			count++;
	
			if(co<actual_no_of_nodes)
			{
				netlist[row][co++]=z;
			}
			else
			{
				co=0;
				row++;
				netlist[row][co++]=z;
			}
			fp1>>buff;

		}
	
	
	float temp_pow=log(nodes)/log(2.0);
	int split1,split2;
	split1=floor(temp_pow/2.0);
	split2=ceil(temp_pow/2.0);
	Graph_inf.No_nodes=nodes;
	Graph_inf.rows=(int)pow(2,split2);
	Graph_inf.colums=(int)pow(2,split1);
	//cout<<Graph_inf.rows<<" here "<<Graph_inf.colums;
	
	fp1.close();	
	int distri=atoi(dist);
	
		if(n2>(float)n1 && distri==1)///
		{
			for(i=actual_no_of_nodes;i<nodes;i++)
			{
				for(j=0;j<nodes;j++)
				{
					if(j>=actual_no_of_nodes)
						netlist[i][j]=0;//netlist[i][j]=MAX;
					else
						netlist[i][j]=0;
				}
			}
			for(i=0;i<actual_no_of_nodes;i++)
			{
				for(j=actual_no_of_nodes;j<nodes;j++)
				{
					netlist[i][j]=0;
				}
			}
			for(i=0;i<nodes;i++)
			{
				netlist[i][i]=0;
			}

		}
	
		else if(n2>(float)n1)
		{
			for(i=actual_no_of_nodes;i<nodes;i++)
			{
				for(j=0;j<nodes;j++)
				{
					if(j>=actual_no_of_nodes)
						netlist[i][j]=MAX;
					else
						netlist[i][j]=0;
				}
			}
			for(i=0;i<actual_no_of_nodes;i++)
			{
				for(j=actual_no_of_nodes;j<nodes;j++)
				{
					netlist[i][j]=0;
				}
			}
			for(i=0;i<nodes;i++)
			{
				netlist[i][i]=0;
			}

		}
	
return netlist;	
		
}


/* Function initialize_add() */
/*Inputs: NO INPUTS */
/*Output: NO OUTPUTS */
		
/* Description : This initializes the Global addresses for a Mesh topology for given number of nodes. It also generates the partition id reference table according to the K-L algorithm. Thus, the K-L algorithm uses this partition ID reference table and does not generate the partition ids at each iteration saving computaional time.  */

void initialize_add()
{
	int nodes=Graph_inf.No_nodes;
	float n4=(log((double)nodes)/log(2.0));
	int n=(int)n4;
	int *row_add,*col_add,*ht_add,i,j,k,**tmp,c,l;
	ht_add =(int*)malloc(nodes*sizeof(int*));
	row_add=(int*)malloc(nodes*sizeof(int*));
	col_add=(int*)malloc(nodes*sizeof(int*));
	
	tmp=(int**)malloc((n-1)*sizeof(int*));
	for(i=0;i<n-1;i++)
		tmp[i]=(int*)malloc(nodes*sizeof(int));
	

	/******* Partition Id generation *******/
	for(i=1;i<n;i++)
	{
		c=0;
		l=0;
		for(j=0;j<(int)pow(2.0,i);j++)
		{
			for(k=0;k<(int)pow(2.0,n-i);k++)
				tmp[i-1][l++]=c;
			c++;
		}

    }
	
	
    /****** Initialization of address ******/
	
	
	for(i=0;i<nodes;i=i+4)
	{
	    row_add[i]=col_add[i]=row_add[i+1]=col_add[i+2]=0;
	    col_add[i+1]=row_add[i+2]=col_add[i+3]=row_add[i+3]=1;
	}
	
	/********* Addressing ***********/
	int count,xyz;
	xyz=1;
	for(i=4;i<nodes;i=i*4)
	{
	    count=0;
	    xyz=xyz*2;
	    for(j=0;((i*count)<nodes);j++)
	    {
	        j=j%4;
	        for(k=0;((k<i) && ((i*count+k)<nodes));k++)
	        {
	            if ((j==2) || (j==3))
                    row_add[i*count+k]+=xyz;
                if ((j==1) || (j==3))
                    col_add[i*count+k]+=xyz;
	        }
	        count++;
	    }
	}
	
	int x,y;
	x = (int)pow(2.0,(n-1)/2);
	if(n%2)
	y=x;
	else
	y=2*x;
	c = 0;

	for(int i = 0; i < nodes; ++i){
		if(i%(nodes/NUM_LAYERS) == 0){
			c++;
		}
		ht_add[i] = c - 1;
	}

	for(int z=0; z<nodes; z++)
	{
		row_add[z]=row_add[z]%x;
		col_add[z]=col_add[z]%y;
		// ht_add[z]=z/(nodes/2);
		cout<<"Z value "<<z<<" ht value "<<ht_add[z]<<"\n";
	}
	//Global address
	address.colum=col_add;
    address.row=row_add;
	address.ht=ht_add;
	print_arr(address.row, nodes);
	print_arr(address.colum, nodes);
	print_arr(address.ht, nodes);

	
    a=tmp;
}



/* Function main() */
/*Inputs: argc => number of command line inputs */
/*		  argv => array of strings ie. 2D array of charactres, it contains the command line inputs. argv[1] holds the name of input graph file, argv[2] holds the initial number of cuts, argv[3] holds if the placement is to be Distributive(1) or Non-Distributive(0) */
/*Output: NO OUTPUTS */
		
/* Description : Topmost function. Controls all the other functions and saves the results*/

/* NOTE : FOR RUNNING THE PROGRAM THE INPUT GRAPHS MUST BE IN THE SAME FOLDER AS THE MAIN AND HEADER FILES. THE WAY OF COMMAND LINE SEQUENCE FOR EXECUTION IS AS FOLLOWS:
./a.out <Graph name> <initial number of cuts> <Distributive(1)/Non-Distributive(0)>
ex.: 	./a.out Graph4.txt 100 0  		*/
int main(int argc, char *argv[])
{
	time_t t;
	time (&t);
	srand(t);
	if (argc < 1)
        {

	  printf("\t\t\t\t\tNO INPUT GRAPH or no. of cuts\n\n");
	  exit(0);
	}
	fstream fptr,fptr1;
	float output[2],**graph;
	graph=Read_graph(argv[1],argv[3]);
	printf("read graph\n");
	char str[30]="final_";
	strcat(str,argv[1]);
	char str1[30]="final_comma_";
	strcat(str1,argv[1]);
	
	fptr.open(str,ios::out | ios::app);
	fptr1.open(str1,ios::out | ios::app);
	fptr1.setf(ios::fixed,ios::floatfield);
	fptr.setf(ios::fixed,ios::floatfield);

	fptr1<<"W,CC,Variance,Peak_Temp\n";

/******power initialization********/
	
	//float *Core_Power;
	//Core_Power=new float[Graph_inf.No_nodes];
	//create_floor_plan();
	//intialize_power_info(Core_Power);
	float result_table [ 2 ];
/*********************************/
	if(argc >= 5){
		NUM_LAYERS = atoi(argv[4]);
	}
	initialize_add();
	int no_cuts = atoi ( argv [ 2 ] );
	
	save_partition = ( int ** ) malloc ( (no_cuts ) * sizeof(int*));
	seed_value = ( int * ) malloc ( (no_cuts ) * sizeof(int));
	
	for ( int j=0; j<no_cuts; j++){
		save_partition [ j ] = ( int * ) malloc ( (Graph_inf.No_nodes ) * sizeof(int)); 
		
	}
	int *core_id = NULL;
	core_id=(int *)malloc(Graph_inf.No_nodes*sizeof(int));
	
	int **partition;
	partition=(int **)malloc(2*sizeof(int*));
	partition[0]=(int *)malloc(Graph_inf.No_nodes/2 * sizeof(int));
	partition[1]=(int *)malloc(Graph_inf.No_nodes/2 * sizeof(int));
	
	for(int j=0;j<Graph_inf.No_nodes;j++)
	{
		core_id[j]=j;
	}
	
	int i;
	int flag=0;
	result_table [ 1 ] = 999999999;
	
		flag=0;
		printf("check\n");
	for( int k = 0; k<NO_OF_RUN; k++)
		{
		flag = 0;
		for ( int j=0; j<no_cuts; j++ ) {
		seed_value [ j ] = rand();
		}
	 	//printf("cuts %d\n", no_cuts);
			perform_KL(argv, graph, output);
			if (result_table [ 1 ] > output [ 1 ] ){
			flag=1;
			}
			if (  flag==1 ){
		  	result_table [ 0 ] = output[0];
		  	result_table [ 1 ] = output[1];
		  	}
		 
			cout<<"\n\n\n\n\nfinished\n\n";
		
	 	
	
		}
		fptr1<<result_table [ 0 ]<<","/*<<result_table [ i ] [ 4 ]<<","<<result_table [ i ] [ 2 ]<<","<<result_table [ i ] [ 3 ]<<","*/<<result_table [ 1 ]<<"\n";
	
	fptr.close();
	fptr1.close();
	return 0;
}


/* Function perform_KL() */
/*Inputs: argv[] => holds the command line inputs */
/*		  netlist => adjacency matrix */
/*		  output => Since there are multiple outputs this array is passed by reference from the controlling function to save results */
/*		  Core_Power => float array holding the power values of the cores */
/*Output: output => float array holding results generated by the algorithm. Passed by reference from controlling function */
		
/* Description : This function is the controlling function of all the partitioning, mapping, cost and thermal (hotspot) functions. Since the results from K-L algorithm depend on initial cut, a number of initial random cuts are generated and then K-L algorithm is called to generate the result for each initial cut. Each of these result is saved and the best result is found out. The best result is returned to the main() function.*/

void perform_KL(char *argv[], float **netlist,float *output)
{

	clock_t t1,t2;
    t1=clock();
    
	int *core_id,i,j,**final_partition_core,actual_no_of_nodes,nodes;
	
	char st1[]="RESULTS628.txt";
	fstream fp;
	fp.open(st1,ios::app | ios::out);
		

	nodes=Graph_inf.No_nodes;
	actual_no_of_nodes=Graph_inf.actual_No_nodes;
	int distri=atoi(argv[3]);
	core_id=(int *)malloc(nodes*sizeof(int));
	for(i=0;i<nodes;i++)
	{
		core_id[i]=i;
	}

	float comm_cost;
	int **partition,init;
	partition=(int **)malloc(2*sizeof(int*));
	partition[0]=(int *)malloc(nodes/2 * sizeof(int));
	partition[1]=(int *)malloc(nodes/2 * sizeof(int));

init=atoi(argv[2]);
/******power variables**************
	float *peak_temp,*var;//,*Core_Power;
	var=new float[init];
	peak_temp= new float[init];
	
/*********************************/

    


	final_partition_core=(int **)malloc(init*sizeof(int*));
	for(i=0;i<init;i++)
		final_partition_core[i]=(int *)malloc(nodes*sizeof(int));
	
float rask;

	for(i=0;i<init;i++)
	{
		
		Init_comm_cut_cost=0;
		Init_powbal=0;
		KL(i,netlist,core_id,nodes,partition);
		
		KL_partition(i,netlist, nodes/2, partition[0],final_partition_core[i]);
		KL_partition(i,netlist, nodes/2, partition[1],final_partition_core[i]);
    		q=0;
rask = cost(final_partition_core[i],netlist,nodes);	
printf("%f",rask);
	for(j = 0; j <nodes; j=j+2)
		{
			fp<<final_partition_core[i][j]+1<<" "<<final_partition_core[i][j+1]+1<<"\t";
		}
		fp<<"\n"<<i<<"\n"<<rask<<"\n\n" ;
	}
	
	float best_cost,temp,*cost_f;
	cost_f= new float[init];
	int best=0;
	
	
	best_cost=999999999;


	for(i=0;i<init;i++)
	{
		cost_f[i]=map_nodes(nodes,final_partition_core[i],netlist);
	}
	
	cout.setf(ios::fixed,ios::floatfield);
	fp.setf(ios::fixed,ios::floatfield);
	
	if(distri==1)
	{
		fp<<"\n"<<"Distributive\t"<<argv[1]<<"\tno. of cores= "<<actual_no_of_nodes<<"\tno. of cuts "<<init;
		fp<<"\n";
		cout<<"\n"<<"Distributive\t"<<argv[1]<<"\tno. of cores= "<<actual_no_of_nodes;
		cout<<"\n";
    }
    else
    {
    	fp<<"\n"<<argv[1]<<"\tno. of cores= "<<actual_no_of_nodes<<"\tno. of cuts "<<init;
		fp<<"\n";
		cout<<"\n"<<argv[1]<<"\tno. of cores= "<<actual_no_of_nodes;
		cout<<"\n";
    }

for(j = 0; j<init;j++)
{
        for(i=0;i<nodes;i=i+2)
        {
     		cout<<final_partition_core[j][i]+1<<" "<<final_partition_core[j][i+1]+1<<"\t";
           	fp<<final_partition_core[j][i]+1<<" "<<final_partition_core[j][i+1]+1<<"\t";
           
        }
cost_f[j] = cost(final_partition_core[j],netlist,nodes);	
temp=cost_f[j];
if(temp<best_cost)
			{
			    best_cost=temp;
			    best=j;
			}
	
fp<<"\n"<<j<<"\n"<<cost_f[j]<<"\n\n";
}

	fp<<"\n\n***** best cost="<<(float)cost_f[best]<<"\nbest ="<<best;

	
	cout<<"\n\n***** best cost=*****"<<cost_f[best]<<"\n";
	output[0]=cost(final_partition_core[best],netlist, nodes);
	cout<<"\n cost in CC : "<<output[0];
	output[1]=cost_f[best];
	
	//cout<<"\nPeak teamperature"<<"\t"<<peak_temp[best]<<"\nVariance\t"<<var[best]<<endl;
		//fp<<"\nPeak teamperature"<<"\t"<<peak_temp[best]<<"\nVariance\t"<<var[best]<<endl;
		//output[1]=var[best];
		//output[2]=peak_temp[best];
		cout<<"\n\n\n*************\n"<<Graph_inf.No_nodes<<"\t"<<Graph_inf.rows<<"\t"<<Graph_inf.colums<<"\n";
	t2=clock();
    float diff =((double) (t2 - t1)) / CLOCKS_PER_SEC;
    fp<<"\tExecution Time: "<<diff<<endl;
    cout<<"\tExecution Time: "<<diff<<endl;

cout<<"\n\n\n\n\nfinished\n\n";
//delete peak_temp;
//delete var;	

free(partition); 
fp.close();


}


/* Function map_nodes() */
/*Inputs: nodes => number of nodes in the graph */
/*		  final_partition_core => array holding the final core sequence generated after partitioning */
/*		  graph => adjacency matrix */
/*		  core_pow => array containing power values of each core */
/*Output: final communication cost of the mapping generated  */
/*		  final_partition_core => It holds the final core sequence generated after different cost and power improvements */
		
/* Description : This function takes the core sequence generated from the partioning algorithm and first applies iterative power improvements (local followed by global) and the iterative cost improvements (local followed by global). Then the final best cost thus found and the corresponding core sequence is returned.  */


float map_nodes(int nodes, int *final_partition_core, float **graph)
{

	int i,j,k;
	int best_partition;
	float temp_cost;

	
	int *temp_final_partition_core;
	temp_final_partition_core=(int*)malloc(nodes*sizeof(int*));
	
	float best_cost[4],Global_best,temp;
 

 /***************** Iterative improvement phase for Communication cost improvements *****************/

/*********************Local phase*********************/   
    
  iterative_improvement(graph,final_partition_core,nodes,1);				//Local swap analysis. '1' stands for local.
    
/******************Global phase *****************************/
 
	 float Final_Global_best_cost=cost(final_partition_core,graph,nodes);
	 
	 for(i=0;i<nodes;i++)
			{
				temp_final_partition_core[i]=final_partition_core[i];
			}
		iterative_improvement(graph,temp_final_partition_core,nodes,0);				
		best_cost[3]=cost(temp_final_partition_core,graph,nodes);
	 	best_cost[1]=best_cost[2]=best_cost[0]=MAX;
	 
	 while(1)
	 {
			iterative_improvement(graph,temp_final_partition_core,nodes,0);
	 	 	temp=cost(temp_final_partition_core,graph,nodes);
	 	 	//cout<<"\nCost: "<<best_cost[3];
	 	 	if(temp<=best_cost[3])
	 	 	{
	 	 		best_cost[0]=best_cost[1];									/* The Global analysis continues till there is no*/ 	 	 		best_cost[1]=best_cost[2];									/* improvement in cost for 4 consecutive passes. */	 	 			best_cost[2]=best_cost[3];									/* The number of comparisions for communication cost*/ 	 	 			best_cost[3]=temp;											/* can be increased if required. 					*/
	
	 	 	}
	 	 	
	 	 	if(temp==best_cost[3] && temp==best_cost[2] &&  temp==best_cost[1] && temp==best_cost[0])
	 	 	{
	 	 		break;
	 	 	}
	 }   
	 
 if(Final_Global_best_cost>best_cost[3])
 {
 	Final_Global_best_cost=best_cost[3];
 	for(i=0;i<nodes;i++)
 			final_partition_core[i]=temp_final_partition_core[i];
 }

Final_Global_best_cost=cost(final_partition_core,graph,nodes);
 
    free(temp_final_partition_core);
    return (Final_Global_best_cost);
}


/* Function iterative_improvement() */
/*Inputs: nodes => number of nodes in the graph */
/*		  final_partition_core => array holding the final core sequence generated after partitioning */
/*		  graph => adjacency matrix */
/*		  local => control varaible which controls the type of pass local or global */
/*Output: final_partition_core => It holds the final core sequence generated after cost improvement */
		
/* Description : This function takes the core sequence generated from the partitioning and performs iterations for improving communication cost of the mapping. The partitions are selected in pairs, at a level of partitioning, and different arrangements are generated in one of the partition keeping the other fixed. The same procedure is repeated for all the levels of partitioning. */

void iterative_improvement(float **graph, int *final_partition_core,int nodes, int local)
{
	float n4=(log((double)nodes)/log(2.0));
	int n=(int)n4;
	int curr_lvl,t,w,min_row,max_row,min_col,max_col,max_ht,min_ht,shift,loop1,best_partition;
    float best_cost[4],avg_row,avg_col,avg_ht,diffrence;
    int i,j,k;
    int level=(int)(log(nodes)/log(2.0))-1;
    float	Global_best;
    
    int **temp_final_partition_core;
	temp_final_partition_core=(int**)malloc(4*sizeof(int*));
	
	temp_final_partition_core[0]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[1]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[2]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[3]=(int*)malloc(nodes*sizeof(int));
	
	int *temp_col,*temp_row,*temp_ht;
	
		
	temp_row=(int*)malloc(nodes*sizeof(int*));
	temp_col=(int*)malloc(nodes*sizeof(int*));
	temp_ht = (int*)malloc(nodes*sizeof(int*));
	for(curr_lvl= level;curr_lvl>0;curr_lvl--)
    {
    
		for(i=0;i<nodes;i++)
		{
			temp_final_partition_core[0][i]=temp_final_partition_core[1][i]=temp_final_partition_core[2][i] 			=temp_final_partition_core[3][i]=final_partition_core[i];
		}
		
	
    	for(i=0;i<(int)pow(2,curr_lvl);i=i+2)
    	{
    		for(j=0;j<nodes;j++)
			{
				temp_row[j]=address.row[j];
				temp_col[j]=address.colum[j];
				temp_ht[j]=address.ht[j];
			}
			
			for(k=0;k<nodes;k++)			
				{
					if(a[curr_lvl-1][k]==i+1)
						break;
				}
			t=(int)pow(2.0,n-curr_lvl);
			
				min_row=address.row[k-t];
				min_col=address.colum[k-t];
				max_row=address.row[k-t];
				max_col=address.colum[k-t];
				max_ht=address.ht[k-t];
				min_ht=address.ht[k-t];
				
				for(w=k-t;w<k;w++)
				{
					if(address.row[w]>max_row)
						max_row=address.row[w];
					
					if(address.colum[w]>max_col)
						max_col=address.colum[w];

					if(address.ht[w]>max_ht)
						max_ht=address.ht[w];
				}
				
				avg_row=(max_row+min_row)/2.0;
				avg_col=(max_col+min_col)/2.0;
				avg_ht=(max_ht+min_ht)/2.0;
			
				best_cost[0]=flip(graph,temp_final_partition_core[0],k,t,nodes,local);
			

			if(curr_lvl == level)
			{
				/**********Flip along lvl axis***********/
				diffrence=(max_ht-min_ht)/2.0;
				shift=ceil(diffrence);
					
					for(w=k-t;w<k;w++)
					{
						if(address.ht[w]>=avg_ht)
							address.ht[w]=address.ht[w]-shift;
						else if(address.ht[w]<avg_ht)
							address.ht[w]=address.ht[w]+shift;
					}
				
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[1][w]=final_partition_core[j];
                            break;
                        }

                    }
                }
                
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				best_cost[1]=flipd(graph,temp_final_partition_core[1],k,t,nodes,local);
				

				
				/***********Flip along Horizontal axis********************************/
				
				diffrence=(max_row-min_row)/2.0;
				shift=ceil(diffrence);
				
				for(w=k-t;w<k;w++)
				{
					if(address.row[w]>=avg_row)
						address.row[w]=address.row[w]-shift;
					else if(address.row[w]<avg_row)
						address.row[w]=address.row[w]+shift;
				}
				
					
					for(w=0;w<nodes;w++)
               		{
                    	for(j=0;j<nodes;j++)
                    	{
                       		 if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        	 {
                            	temp_final_partition_core[2][w]=final_partition_core[j];
                            	break;
                        	 }
                    	}
                	}
                	
                	for(j=0;j<nodes;j++)
					{
						address.row[j]=temp_row[j];
						address.colum[j]=temp_col[j];
						address.ht[j]=temp_ht[j];
					}
					
				best_cost[2]=flipd(graph,temp_final_partition_core[2],k,t,nodes,local);	
        
		/**********Flip along lvl axis***********/	
		
				diffrence=(max_ht-min_ht)/2.0;
				shift=ceil(diffrence);
					
					for(w=k-t;w<k;w++)
					{
						if(address.ht[w]>=avg_ht)
							address.ht[w]=address.ht[w]-shift;
						else if(address.ht[w]<avg_ht)
							address.ht[w]=address.ht[w]+shift;
					}
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[3][w]=final_partition_core[j];
                            break;
                        }

                    }
                }
                
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
			best_cost[3]=flip(graph,temp_final_partition_core[3],k,t,nodes,local);	
			
	}
	else
	{
			/**********Flip along horizontal axis***********/
				diffrence=(max_row-min_row)/2.0;
				shift=ceil(diffrence);
					
					for(w=k-t;w<k;w++)
					{
						if(address.row[w]>=avg_row)
							address.row[w]=address.row[w]-shift;
						else if(address.row[w]<avg_row)
							address.row[w]=address.row[w]+shift;
					}
				
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[1][w]=final_partition_core[j];
                            break;
                        }

                    }
                }
                
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				best_cost[1]=flipd(graph,temp_final_partition_core[1],k,t,nodes,local);
				

				
				/***********Flip along Vertical axis********************************/
				
				diffrence=(max_col-min_col)/2.0;
				shift=ceil(diffrence);
				
				for(w=k-t;w<k;w++)
				{
					if(address.colum[w]>=avg_col)
						address.colum[w]=address.colum[w]-shift;
					else if(address.colum[w]<avg_col)
						address.colum[w]=address.colum[w]+shift;
				}
				
					
					for(w=0;w<nodes;w++)
               		{
                    	for(j=0;j<nodes;j++)
                    	{
                       		 if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        	 {
                            	temp_final_partition_core[2][w]=final_partition_core[j];
                            	break;
                        	 }
                    	}
                	}
                	
                	for(j=0;j<nodes;j++)
					{
						address.row[j]=temp_row[j];
						address.colum[j]=temp_col[j];
						address.ht[j]=temp_ht[j];
					}
					
				best_cost[2]=flip(graph,temp_final_partition_core[2],k,t,nodes,local);	
        
		/**********Flip along horizontal axis***********/	
		
				diffrence=(max_row-min_row)/2.0;
				shift=ceil(diffrence);
					
					for(w=k-t;w<k;w++)
					{
						if(address.row[w]>=avg_row)
							address.row[w]=address.row[w]-shift;
						else if(address.row[w]<avg_row)
							address.row[w]=address.row[w]+shift;
					}
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[3][w]=final_partition_core[j];
                            break;
                        }

                    }
                }
                
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
			best_cost[3]=flip(graph,temp_final_partition_core[3],k,t,nodes,local);	
				
	}
				/*******Checked all combinations*******/
			
			Global_best=best_cost[0];	
			best_partition=0;
			for(j=0;j<4;j++)
			{
				if(Global_best>best_cost[j])
				{
					Global_best=best_cost[j];
					best_partition=j;
				}
			}

			
			for(j=0;j<nodes;j++)
                {
                	final_partition_core[j]=temp_final_partition_core[best_partition][j];
                }

		}
    }
   
   
    free(temp_col);
    free(temp_row);
   	free(temp_final_partition_core); 
    //cout<<"\nIterative Improvement Phase Complete\n";
   
}


/* Function flip() */
/*Inputs: nodes => number of nodes in the graph */
/*		  final_partition_core => array holding the final core sequence generated after partitioning */
/*		  G => adjacency matrix */
/*		  local => control varaible which controls the type of pass local or global */
/*		  k => start index of the partition under consideration*/
/*		  t => size of the partition under consideration*/
/*Output: final_partition_core => It holds the final core sequence generated after cost improvement */
		
/* Description : This function takes a partition and rearranges the core within it and checks for cost improvements. The arrangement with the lowest cost is returned.*/

float flip(float **G,int *final_partition_core,int k,int t,int nodes,int local)
{

	float cost_arr[4],best_cost,avg_row,avg_col,avg_ht,diffrence;
	if(local==1)
		cost_arr[0]=best_cost=cost_local(final_partition_core,G,nodes,k-t,k+t);
	else if(local==0)
		cost_arr[0]=best_cost=cost(final_partition_core,G,nodes);
		    	
	int **temp_final_partition_core,best=0,w,min_row,max_row,min_col,max_col,min_ht,max_ht,shift,j,*temp_col,*temp_row,*temp_ht;
	temp_final_partition_core=(int**)malloc(4*sizeof(int*));
	
	temp_final_partition_core[0]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[1]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[2]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[3]=(int*)malloc(nodes*sizeof(int));
	
	temp_row=(int*)malloc(nodes*sizeof(int*));
	temp_col=(int*)malloc(nodes*sizeof(int*));
	temp_ht=(int*)malloc(nodes*sizeof(int*));
	    	for(j=0;j<nodes;j++)
			{
				temp_row[j]=address.row[j];
				temp_col[j]=address.colum[j];
				temp_ht[j]=address.ht[j];
				temp_final_partition_core[0][j]=final_partition_core[j];
			}
			
			
				min_row=address.row[k];
				min_col=address.colum[k];
				max_row=address.row[k];
				max_col=address.colum[k];
				min_ht=address.ht[k];
				max_ht=address.ht[k];
				for(w=k;w<k+t;w++)
				{
					if(address.row[w]>max_row)
						max_row=address.row[w];
						
					if(address.colum[w]>max_col)
						max_col=address.colum[w];

					if(address.ht[w]>max_ht)
						max_ht=address.ht[w];
				}
				avg_col=(max_col+min_col)/2.0;
				avg_row=(max_row+min_row)/2.0;
				avg_ht=(max_ht+min_ht)/2.0;
				
				/********Flip along horizontal axis*********/
				
				
				diffrence=(max_row-min_row)/2.0;
				shift=ceil(diffrence);
       			for(w=k;w<k+t;w++)
				{
					if(address.row[w]>=avg_row)
						address.row[w]=address.row[w]-shift;
					else if(address.row[w]<avg_row)
						address.row[w]=address.row[w]+shift;
				}
				
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[1][w]=final_partition_core[j];

                            //temp_final_partition_core[w+1]=final_partition_core[j+1];
                            break;
                        }

                    }
                }
                	
			
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				if(local==1)
					cost_arr[1]=cost_local(temp_final_partition_core[1],G,nodes,k-t,k+t);
				else
					cost_arr[1]=cost(temp_final_partition_core[1],G,nodes);
				
				
		/********Flip along vertical axis*********/
;			
				
				
				diffrence=(max_col-min_col)/2.0;
				shift=ceil(diffrence);
                   
				for(w=k;w<k+t;w++)
				{
					if(address.colum[w]>=avg_col)
						address.colum[w]=address.colum[w]-shift;
					else if(address.colum[w]<avg_col)
						address.colum[w]=address.colum[w]+shift;
				}
																				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[2][w]=final_partition_core[j];
                            //temp_final_partition_core[w+1]=final_partition_core[j+1];
                            break;
                        }

                    }
                }
                
       
				
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				if(local==1)
					cost_arr[2]=cost_local(temp_final_partition_core[2],G,nodes,k-t,k+t);
				else
					cost_arr[2]=cost(temp_final_partition_core[2],G,nodes);
	
	
	/*******************Flip along horizontal axis*********************/
	
	
				diffrence=(max_row-min_row)/2.0;
				shift=ceil(diffrence);
       
				for(w=k;w<k+t;w++)
				{
					if(address.row[w]>=avg_row)
						address.row[w]=address.row[w]-shift;
					else if(address.row[w]<avg_row)
						address.row[w]=address.row[w]+shift;
				}
				
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[3][w]=final_partition_core[j];
                            //temp_final_partition_core[w+1]=final_partition_core[j+1];
                            break;
                        }

                    }
                }
                
       
				
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				if(local==1)
					cost_arr[3]=cost_local(temp_final_partition_core[3],G,nodes,k-t,k+t);
				else
					cost_arr[3]=cost(temp_final_partition_core[3],G,nodes);
					
				
					
	/*************Finding the best partitioning****************************/
					
				for(j=1;j<4;j++)
				{
					if(cost_arr[j]<best_cost)
					{
						best_cost=cost_arr[j];
						best=j;
					}
				}
					
	/**************** Finalizing results ***********************/
				for(j=0;j<nodes;j++)
                {
                	final_partition_core[j]=temp_final_partition_core[best][j];
                }
                   
            
            free(temp_col);
    	free(temp_row);
  	 free(temp_final_partition_core);    
                
                
	return best_cost;				
}

float flipd(float **G,int *final_partition_core,int k,int t,int nodes,int local)
{

	float cost_arr[4],best_cost,avg_row,avg_col,avg_ht,diffrence;
	if(local==1)
		cost_arr[0]=best_cost=cost_local(final_partition_core,G,nodes,k-t,k+t);
	else if(local==0)
		cost_arr[0]=best_cost=cost(final_partition_core,G,nodes);
		    	
	int **temp_final_partition_core,best=0,w,min_row,max_row,min_col,max_col,min_ht,max_ht,shift,j,*temp_col,*temp_row,*temp_ht;
	temp_final_partition_core=(int**)malloc(4*sizeof(int*));
	
	temp_final_partition_core[0]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[1]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[2]=(int*)malloc(nodes*sizeof(int));
	temp_final_partition_core[3]=(int*)malloc(nodes*sizeof(int));
	
	temp_row=(int*)malloc(nodes*sizeof(int*));
	temp_col=(int*)malloc(nodes*sizeof(int*));
	temp_ht=(int*)malloc(nodes*sizeof(int*));
	    	for(j=0;j<nodes;j++)
			{
				temp_row[j]=address.row[j];
				temp_col[j]=address.colum[j];
				temp_ht[j]=address.ht[j];
				temp_final_partition_core[0][j]=final_partition_core[j];
			}
			
			
				min_row=address.row[k];
				min_col=address.colum[k];
				max_row=address.row[k];
				max_col=address.colum[k];
				min_ht=address.ht[k];
				max_ht=address.ht[k];
				for(w=k;w<k+t;w++)
				{
					if(address.row[w]>max_row)
						max_row=address.row[w];
						
					if(address.colum[w]>max_col)
						max_col=address.colum[w];

					if(address.ht[w]>max_ht)
						max_ht=address.ht[w];
				}
				avg_col=(max_col+min_col)/2.0;
				avg_row=(max_row+min_row)/2.0;
				avg_ht=(max_ht+min_ht)/2.0;
				
				/********Flip along lvl axis*********/
				
				
				diffrence=(max_ht-min_ht)/2.0;
				shift=ceil(diffrence);
       			for(w=k;w<k+t;w++)
				{
					if(address.ht[w]>=avg_ht)
						address.ht[w]=address.ht[w]-shift;
					else if(address.ht[w]<avg_ht)
						address.ht[w]=address.ht[w]+shift;
				}
				
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[1][w]=final_partition_core[j];

                            //temp_final_partition_core[w+1]=final_partition_core[j+1];
                            break;
                        }

                    }
                }
                	
			
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				if(local==1)
					cost_arr[1]=cost_local(temp_final_partition_core[1],G,nodes,k-t,k+t);
				else
					cost_arr[1]=cost(temp_final_partition_core[1],G,nodes);
				
				
		/********Flip along horizontal axis*********/
;			
				
				
				diffrence=(max_row-min_row)/2.0;
				shift=ceil(diffrence);
                   
				for(w=k;w<k+t;w++)
				{
					if(address.row[w]>=avg_row)
						address.row[w]=address.row[w]-shift;
					else if(address.row[w]<avg_row)
						address.row[w]=address.row[w]+shift;
				}
																				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[2][w]=final_partition_core[j];
                            //temp_final_partition_core[w+1]=final_partition_core[j+1];
                            break;
                        }

                    }
                }
                
       
				
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				if(local==1)
					cost_arr[2]=cost_local(temp_final_partition_core[2],G,nodes,k-t,k+t);
				else
					cost_arr[2]=cost(temp_final_partition_core[2],G,nodes);
	
	
	/*******************Flip along lvl axis*********************/
	
	
				diffrence=(max_ht-min_ht)/2.0;
				shift=ceil(diffrence);
       
				for(w=k;w<k+t;w++)
				{
					if(address.ht[w]>=avg_ht)
						address.ht[w]=address.ht[w]-shift;
					else if(address.ht[w]<avg_ht)
						address.ht[w]=address.ht[w]+shift;
				}
				
				
				for(w=0;w<nodes;w++)
                {
                    for(j=0;j<nodes;j++)
                    {
                        if(temp_row[w]==address.row[j] && temp_col[w]==address.colum[j] && temp_ht[w]==address.ht[j])
                        {
                            temp_final_partition_core[3][w]=final_partition_core[j];
                            //temp_final_partition_core[w+1]=final_partition_core[j+1];
                            break;
                        }

                    }
                }
                
       
				
				for(j=0;j<nodes;j++)
				{
					address.row[j]=temp_row[j];
					address.colum[j]=temp_col[j];
					address.ht[j]=temp_ht[j];
				}
				
				if(local==1)
					cost_arr[3]=cost_local(temp_final_partition_core[3],G,nodes,k-t,k+t);
				else
					cost_arr[3]=cost(temp_final_partition_core[3],G,nodes);
					
				
					
	/*************Finding the best partitioning****************************/
					
				for(j=1;j<4;j++)
				{
					if(cost_arr[j]<best_cost)
					{
						best_cost=cost_arr[j];
						best=j;
					}
				}
					
	/**************** Finalizing results ***********************/
				for(j=0;j<nodes;j++)
                {
                	final_partition_core[j]=temp_final_partition_core[best][j];
                }
                   
            
            free(temp_col);
    	free(temp_row);
  	 free(temp_final_partition_core);    
                
                
	return best_cost;				
}



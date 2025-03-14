#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random>
#include <unordered_map>
#include <set>
#include "global.h"

using namespace std;


//Declare global functions
float** allocateMatrix(int rows, int cols);
void deallocateMatrix(float** matrix, int rows);


void writeToFile(unordered_map<int, float>& values, int type, int iter){
	int s = values.size();
	string filename = "";
	if(type == 1){
		filename = "iteration_"+to_string(iter)+"_pop_";
	}else{
		filename = "lifetime_values";
	}
	filename = filename+to_string(s)+".csv";
	ofstream file(filename);

	if(!file.is_open()){
		cerr<< "Error opening file: "<<filename << endl;
	}

	file << "Worker,value\n";
	for(int i=0;i<values.size();++i){
		float v = values[i];
		file<< i << "," << v <<"\n";
	}
}

void drawToFile(float** &network, int population, int iter){
	string filename = "graph"+to_string(iter)+".txt";
	ofstream file(filename);

	if(!file.is_open()){
		cerr<< "Error opening file: "<<filename << endl;
	}

	for(int i=0;i<population;++i){
		for(int j=0;j<population;++j){
			file<< network[i][j] << " ";
		}
		file<< "\n";
	}
}
set<int> getNeighbors(float** &network, int population, int node){
	set<int> neighbors;
	for(unsigned int i=1;i<population;i++){
		if(network[node][i] != 0 && node != i){
			neighbors.insert(i);
		}
	}

	return neighbors;
}
/*
 *Return a stochastic matrix for propagating through the network. 
NOTE: This puts zeros on the first row/col to avoid propagating back to the queen
 * */
float** &getStochasticMatrix (float** network, unsigned int population, unordered_map<int, int> &degrees, float** &stochastic_matrix){
	for(int i=0;i<population;++i){
		float p = 0.9;
		for(int j=0;j<population;++j){
			p = (j==0)?1:0.9; //Share everything from the queen.
			int x = degrees[j];
			if(x==0){
				p = 0;//If no neighbors, hold on to qmp
			}
			if( network[i][j] != 0 && i>0 ){
				float entry = (i==j) ? 1-p:(p/x);//*network[i][j]; 
				stochastic_matrix[i][j] = entry;
				
			}else{
				stochastic_matrix[i][j] = 0;
				stochastic_matrix[j][i] = 0;
			}
		}

	}

	return stochastic_matrix;

}

void changeRetinue(float** network, unsigned int population, unordered_map<int, float> &values, float threshold, int* &ret, int ret_size, random_device &rd, unordered_map<int, int> &deg, unordered_map<int, bool> &isRet){
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(1, population-1);
	int newMember = dis(gen);
	int counter = 0;
	while((values[newMember] < threshold && counter<population) || isRet[newMember]){
		newMember = dis(gen);
		counter += 1;
	}
	uniform_int_distribution<int> dis2(0, ret_size-1);
	int oldMemberIndex = dis2(gen);
	int duck = dis(gen);

	for(int k=0;k<ret_size;++k){
		if(ret[k] != ret[oldMemberIndex]){
			network[ret[k]][ret[oldMemberIndex]] = 0;
			network[ret[oldMemberIndex]][ret[k]] = 0;
			deg[ret[oldMemberIndex]] -= 1;
		}

		uniform_real_distribution<float> str(0.090,0.099);
		float connection = str(gen)+0.9;
		network[ret[k]][newMember] = connection;
		network[newMember][ret[k]] = connection;
		deg[newMember] += 1;
	}
	network[0][newMember] = network[0][ret[oldMemberIndex]];
	network[newMember][0] = network[ret[oldMemberIndex]][0];
	network[0][ret[oldMemberIndex]] = 0;
	network[ret[oldMemberIndex]][0] = 0;
	isRet[ret[oldMemberIndex]] = false;
	isRet[newMember] = true;
	ret[oldMemberIndex] = newMember;

}




/*
 * **********MATRIX OPERATIONS***************************
 *
 *Matrix multiplication. Check that dimensions are compatible before calling
 * */
float** MATMUL(float** A, float** B,  int m, int n, int p){
	float** C = allocateMatrix(m, p);
	for(int i=0;i<m;++i){
		for(int j=0;j<p;++j){
			C[i][j] = 0;
			for(int k=0;k<n;++k){
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}

float** allocateMatrix(int rows, int cols){
	float** matrix = new float*[rows];
	for(int i=0;i<rows;++i){
		matrix[i] = new float[cols];
	}
	return matrix;
}

void deallocateMatrix(float** matrix, int rows){
	for(int i=0;i<rows;++i){
		delete[] matrix[i];
	}
	delete[] matrix;
}

float getConnectionStrength(random_device &rd){
	mt19937 gen(rd());
	uniform_real_distribution<float> str(0.090,0.099);
	return (str(gen)+0.9);
}


void connectGraph(float** &network, int population, int maxDegree, unordered_map<int, int> &degrees, int* &retinue, int retSize, int cliqueSize, random_device &rd, unordered_map<int, bool> &isRet){

	//Clear Graph
	for(int i=0;i<population;i++){
		for(int j=0;j<population;j++){
			network[i][j] = 0;
			if(i==j){
				network[i][j] = 1;
			}
		}
		degrees[i] = 0;
	}


	//Connect Queen to Retinue
	for(int k=0;k<retSize;k++){
		int w = retinue[k];
		float connection = getConnectionStrength(rd);
		network[0][w] = network[w][0] = connection;
	//	degrees[w] += 1;
		degrees[0] += 1;
	//	qmp_values[w] = produced_qmp/retinue_size;
	}

	//Connect Retinue to each other
	for(int i=0;i<retSize;i++){
		for(int j=i+1;j<retSize;j++){
			float connection = getConnectionStrength(rd);
			network[retinue[j]][retinue[i]] = connection;
			network[retinue[i]][retinue[j]] = connection;
			degrees[retinue[i]] += 1;
			degrees[retinue[j]] += 1;
			}
	}

	//Connect rest of graph
	int start = 1;
	int end = start+cliqueSize;
	float p = 0.5f;
	mt19937 gen(rd());
	while(start<population){
		
		uniform_int_distribution<int> pop(start, end-1);
		for(int w=start;w<end;++w){		
			int max = (isRet[w])?(maxDegree*2):maxDegree;
			int deg = (max - degrees[w])*p;
			for(int j=0;j<deg;j++){
				int n = pop(gen);
				if(network[w][n] == 0){				
					int neighborMax = (isRet[n])?maxDegree*2:maxDegree;
					if(degrees[n] < neighborMax){
						float connection = getConnectionStrength(rd);
						network[w][n] = connection;
						network[n][w] = connection;
						degrees[w] = degrees[w]+1;
						degrees[n] = degrees[n]+1;
					}
				}
			}

			max = (isRet[w])?(maxDegree*2):maxDegree;
			deg = (max - degrees[w])*p/2;
			uniform_int_distribution<int> dis2(0, population);
			for(int r=0;r<deg;++r){
				int neighbor = dis2(gen);
				int counter = 0;
				while(network[neighbor][w] != 0 && counter<100){
					neighbor = dis2(gen);
					counter += 1;	
				}
				if(network[neighbor][w] == 0){
					network[w][neighbor] = network[neighbor][w] = getConnectionStrength(rd);
				degrees[w] += 1;
				degrees[neighbor] += 1;
				}
			}
		}
		start += cliqueSize;
		end = ((start+cliqueSize)<population)?start+cliqueSize:population;
	}
}

//argv[??, max_network_size, population, iterations, prop_cycles, max_degree]
int main(int argc, char* argv[]){
	
	const unsigned int initial_population = atoi(argv[2]);	
	unsigned int max_degree = atoi(argv[5]);
	unsigned int retinue_size = 12;
	float produced_qmp = 1000;
	const float THRESH = 0.5;
	unordered_map<int, float> qmp_values;
	unordered_map<int, int> degrees;
	unordered_map<int, bool> isRet;
	unordered_map<int, int> quorum;
	unordered_map<int, float> lifetime_values;
	unsigned int population = initial_population;
        const size_t network_rows = atoi(argv[1]);
	const size_t network_cols = network_rows;
	const int max_iterations = atoi(argv[3]);
	int prop_cycles = atoi(argv[4]);
	/*Create Large Blank Network */
	float** network = allocateMatrix(network_rows, network_cols);
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<> pop(1, population-1);
	
	unsigned int cohort_0 = 1;
	unsigned int cohort_n = population-1;
	
	int cliqueSize = population/10;
	
	uniform_int_distribution<> ret(cohort_0, cohort_n);
	uniform_real_distribution<float> str(0.090,0.099);
	int* retinue = new int[retinue_size];


		
	float** stochastic_matrix = allocateMatrix(population, population);
	float** initial_values = allocateMatrix(population, 1);
	for(unsigned int i=0;i<population;i++){
		isRet[i] = false;
	}
	/*Pick initial retinue */
	for(unsigned int k=0;k<retinue_size;k++){
		int w = ret(gen);
		while(network[0][w] != 0){
			w = ret(gen);
		}
		retinue[k] = w;
		isRet[w] = true;
	//	qmp_values[w] = produced_qmp/retinue_size;
	}
	int simIter = 0;
	int quorum_size = 0;
	while(simIter<max_iterations && quorum_size<(population/10)){
	/*Clear Network*/
		for(unsigned int i=0;i<population;i++){
			qmp_values[i] = 0;
			lifetime_values[i] = 0;
		}

	/*Connect retinue to each other*/

	/*Connect rest of network*/


	 	connectGraph(network, population, max_degree, degrees, retinue, retinue_size, cliqueSize, rd, isRet);

//		drawToFile(network, population, simIter);
		/*Propagate*/
		int step = 0;
		float threshold = 0;
	//	stochastic_matrix = getStochasticMatrix(network, population, degrees, stochastic_matrix);
		while(step<prop_cycles){
			stochastic_matrix = getStochasticMatrix(network, population, degrees, stochastic_matrix);
			qmp_values[0] = produced_qmp/prop_cycles;
			for(int i=0;i<population;++i){
				initial_values[i][0] = qmp_values[i];
			}
	/**************Propagation through network****************************/	
			initial_values = MATMUL(stochastic_matrix, initial_values, population, population, 1);
	
	/********************************************************************/
			for(int i=0;i<population;++i){
				qmp_values[i] = initial_values[i][0];
				lifetime_values[i] += qmp_values[i];
				unordered_map<int,int>::iterator it = quorum.find(i);
				if(qmp_values[i] < THRESH){
					if(it != quorum.end()){
						quorum[i] += 1;
					}else{
						quorum[i] = 1;
					}
				}else{
					if(it != quorum.end()){
						if(quorum[i] < prop_cycles){
							quorum.erase(i);
						}
					}
				}

			}
	//	deallocateMatrix(stochastic_matrix, population);
			if(step>0 && step%3 == 0){
				changeRetinue(network, population, qmp_values, threshold,  retinue, retinue_size, rd, degrees, isRet);
				threshold = (threshold<THRESH)?0.02*step:THRESH;
			//step%1 will not converge quickly	
	 			if(step%2 == 0){
					connectGraph(network, population, max_degree, degrees, retinue, retinue_size, cliqueSize, rd, isRet);
				}
			}
			step += 1;
	}

		cout << "******Iteration: " << simIter << "***********" <<endl;
	/*Print Netwwork*/
	/*
	for(unsigned int i=0;i<population;i++){
		for(unsigned int j=0;j<population;j++){
			cout << network[i][j] << " ";
		}
		cout << endl;
	}

	cout <<endl;
	*/
	/*Print qmp values*/
		int total_in_system = 0;
		quorum_size = 0;
		cout << "Worker  " <<endl;
		for(unsigned int k=0;k<population;k++){
			cout << k << ": " << qmp_values[k];
			if(quorum.find(k) != quorum.end()){
				int q = quorum[k]/prop_cycles;
				if(q>0){
					cout << " q --> "<< q;
					quorum_size += 1;
				}
			}
			cout << endl;
			total_in_system += qmp_values[k];
		}
		cout << endl;
		cout<< "Retinue" << endl;
		for(int i=0;i<retinue_size;++i){
			cout<< retinue[i] <<" : "<< qmp_values[retinue[i]]<<endl; 
		}
		cout<< endl;
		cout<<"total qmp in system: "<<total_in_system<<endl;
		cout<< "quorum size: "<<quorum_size<<endl;
		
		cout<<endl;
		int half = max_iterations/2;
		if(simIter%half == 0){
			writeToFile(qmp_values, 1, simIter);
		}
		simIter += 1;
	}


	writeToFile(qmp_values, 1, simIter);
	writeToFile(lifetime_values, 2, 0);
//	drawToFile(network, population, max_iterations);
	/*Free network memory*/
	deallocateMatrix(initial_values, population);
	deallocateMatrix(stochastic_matrix, population);
	deallocateMatrix(network, network_rows);
	delete[] retinue;
	return 0;
}

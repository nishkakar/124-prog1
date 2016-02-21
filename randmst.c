#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

float rand_float() {
    return (float) rand() / (float) RAND_MAX;
}

float** init_graph(int numpoints) {
	// allocate graph
	float** graph = (float**) malloc(numpoints * sizeof(float*));
	for (int i = 0; i < numpoints; ++i) {
		graph[i] = (float*) malloc(numpoints * sizeof(float));
		memset(graph[i], 0, numpoints * sizeof(float));
	}

	return graph;
}

int find(int vertex, int* parents) {
	if (vertex != parents[vertex]) {
		parents[vertex] = find(parents[vertex], parents);
	}

	return parents[vertex];
}

void link(int component1, int component2, int* parents, int* ranks) {
	if (ranks[component1] > ranks[component2]) {
        parents[component2] = component1;
	}
	else if (ranks[component1] < ranks[component2]) {
		parents[component1] = component2; 
	}
	else if (ranks[component1] == ranks[component2]) {
		++ranks[component1];
		parents[component2] = component1;		
	}
}

void ds_union(int vertex1, int vertex2, int* parents, int* ranks) {
	link(find(vertex1, parents), find(vertex2, parents), parents, ranks);
}

int CompareArrays(const void* arr1, const void* arr2) {
    // convert to correct type
    const float* one = (const float*) arr1;
    const float* two = (const float*) arr2;

    // only compare first element of array
    if (one[0] < two[0]) return -1;
    if (one[0] > two[0]) return +1;

    // edge weights are the same
    return 0;
}

float** find_mst(float** graph, int numpoints) {
    float** mst = init_graph(numpoints);

	int numedges = ((numpoints) * (numpoints - 1))/2;

	float** sorted_edges = malloc(numedges * sizeof(float*) + (numedges * 3 * sizeof(float)));
	bool contiguous = false;
	if (sorted_edges != NULL) {
	    float* pos = (float*) (sorted_edges + numedges);
	    for (int i = 0; i < numedges; ++i) {
	        sorted_edges[i] = pos + i * 3;
	    }
	    contiguous = true;
	}
	else {
		sorted_edges = malloc(numedges * sizeof(float*));
		for (int i = 0; i < numedges; ++i) {
			sorted_edges[i] = malloc(3 * sizeof(float));
		}
	}

	// fill sorted edges with unsorted edges and edge weights from graph
	int ctr = 0;
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < i; ++j) {
			sorted_edges[ctr][0] = graph[i][j];
			sorted_edges[ctr][1] = i;
			sorted_edges[ctr][2] = j;
			++ctr;
		}
    }

	// sort edges based on edge weight
	qsort(*sorted_edges, numedges, sizeof(float[3 ]), CompareArrays);

    // make sets
    int* parents = (int*) malloc(numpoints * sizeof(int));
    int* ranks = (int*) malloc(numpoints * sizeof(int));
    for (int i = 0; i < numpoints; ++i) {
        parents[i] = i;
        ranks[i] = 1;
    }

    for (int i = 0; i < numedges; ++i) {
        int u = sorted_edges[i][1];
        int v = sorted_edges[i][2];
        if (find(u, parents) != find(v, parents)) {
            mst[u][v] = sorted_edges[i][0];
            mst[v][u] = sorted_edges[i][0];
            ds_union(u, v, parents, ranks);
        }
    }

    // free all malloced memory
    if (!contiguous) {
    	for (int i = 0; i < numedges; ++i) {
    		free(sorted_edges[i]);
    	}
    }
    free(sorted_edges);
    free(parents);
    free(ranks);

	return mst;
}

float** generate_complete_graph(int numpoints) {
	float** graph = init_graph(numpoints);

	// fill graph with random numbers
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < i; ++j) {
			graph[i][j] = rand_float();
		}
	}

	return graph;
}

float distance(float* point1, float* point2, int d) {
	double sum = 0;
	for (int i = 0; i < d; ++i) {
		sum += pow((double) point1[i] - point2[i], (double) 2);
	}
	return pow(sum, 0.5);
}

float** generate_euclidean_graph(int numpoints, int dim) {
	float** graph = init_graph(numpoints);

	// temporary storage for coordinates
	float coordinates[numpoints][dim];
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < dim; ++j) {
			coordinates[i][j] = rand_float();
		}
	}

	// fill graph with random euclidean distances
	// by randomly generating coordinates
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < i; ++j) {
			graph[i][j] = distance(coordinates[i], coordinates[j], dim);
		}
	}

	return graph;
}

void print_graph(float** graph, int numpoints) {
	// prints out adjacency matrix
    for (int i = 0; i < numpoints; ++i) {
    	for (int j = 0; j < numpoints; ++j) {
    		printf("%6.4f ", graph[i][j]);
    	}
    	printf("\n");
    }
}

void free_graph(float** graph, int numpoints) {
	for (int i = 0; i < numpoints; ++i) {
		free(graph[i]);
	}
	free(graph);
}

float find_mst_weight(float** graph, int numpoints) {
    float total_weight = 0.0;
    for (int i = 0; i < numpoints; ++i) {
        for (int j = 0; j < i; ++j) {
            total_weight += graph[i][j];
        }
    }
    return total_weight;
}

int main(int argc, char *argv[]) {
    // // seed rand() w/ current time
    time_t t;
    srand((unsigned) time(&t));

    int numpoints = atoi(argv[2]);
    int iterations = atoi(argv[3]);
    int dimension = atoi(argv[4]);


    printf("%s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4]);

    float average_weight, average_time;

    if (dimension == 0) {
        float complete_mst_weights = 0.0;
        float total_time = 0.0;
        for (int i = 0; i < iterations; ++i) {
            float** complete_graph = generate_complete_graph(numpoints);

            clock_t start = clock(); 
            float** complete_mst = find_mst(complete_graph, numpoints);
            
            total_time += (float) (clock() - start) / CLOCKS_PER_SEC;

            complete_mst_weights += find_mst_weight(complete_mst, numpoints);

            free_graph(complete_graph, numpoints);
            free_graph(complete_mst, numpoints);
        }
        average_weight = complete_mst_weights / (float) iterations;
        average_time = total_time / (float) iterations;
    }
    else {
        float euc_mst_weights = 0.0;
        float total_time = 0.0;
        for (int i = 0; i < iterations; ++i) {
            float** euc_graph = generate_euclidean_graph(numpoints, dimension);
            
            clock_t start = clock(); 
            float** euc_mst = find_mst(euc_graph, numpoints);
            
            total_time += (float) (clock() - start) / CLOCKS_PER_SEC;

            euc_mst_weights += find_mst_weight(euc_mst, numpoints);

            free_graph(euc_graph, numpoints);
            free_graph(euc_mst, numpoints);
        }
        average_weight = euc_mst_weights/ (float) iterations;
        average_time = total_time / (float) iterations;
    }

    printf("Average weight: %f\n", average_weight);
    printf("Average time: %f\n", average_time);

    printf("==========\n");
    return 0;
}

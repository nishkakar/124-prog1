#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// returns random float between 0.0 and 1.0
float rand_float() {
    return (float) rand() / (float) RAND_MAX;
}

/* initializes a complete graph with numpoints # of vertices;
   we implemented the graph in a form of an adjacency matrix */
float** init_graph(int numpoints) {
	float** graph = (float**) malloc(numpoints * sizeof(float*));
	for (int i = 0; i < numpoints; ++i) {
		graph[i] = (float*) malloc(numpoints * sizeof(float));
		memset(graph[i], 0, numpoints * sizeof(float));
	}
	return graph;
}

// frees all memory allocated for adjacency matrix graph representation
void free_graph(float** graph, int numpoints) {
	for (int i = 0; i < numpoints; ++i) {
		free(graph[i]);
	}
	free(graph);
}

/* recursively finds the parent node of a given vertex in order to 
   return which component the vertex belongs to */
int find(int vertex, int* parents) {
	if (vertex != parents[vertex]) {
		parents[vertex] = find(parents[vertex], parents);
	}

	return parents[vertex];
}

/* links two components based on rank; assigns the higher ranked component
   to be the parent of the lower ranked component to ensure that the tree
   is balanced; if ranks are equal, defaults to component1 being assigned as
   component2's parent and increments the rank of component1. */
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

/* given two vertexes, finds their components and links them */
void ds_union(int vertex1, int vertex2, int* parents, int* ranks) {
	link(find(vertex1, parents), find(vertex2, parents), parents, ranks);
}

// comparator we pass into qsort in order to sort edges in our adjacency matrix by weight
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

// implemented Kruskal's algorithm to find the minimum spanning tree in a given graph
float** find_mst(float** graph, int numpoints, unsigned int* numedges_ptr) {
	unsigned int numedges = *numedges_ptr;

    // allocate memory for edges sorted by weight; store the weight and the two vertexes of the edge
	float** sorted_edges = malloc(numedges * sizeof(float*) + (numedges * 3 * sizeof(float)));
    
    /* if we were able to allocate all necessary memory at once, assigns memory;
       else makes separate malloc calls edge by edge */
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

	// fill sorted_edges with unsorted edges and edge weights from graph
	int ctr = 0;
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < i; ++j) {
			if (graph[i][j] != 0.0) {
				sorted_edges[ctr][0] = graph[i][j];
				sorted_edges[ctr][1] = i;
				sorted_edges[ctr][2] = j;
				++ctr;
			}
		}
    }

    // since we don't need this graph we free the memory to allow space for the mst graph
    free_graph(graph, numpoints);
    float** mst = init_graph(numpoints);

	// sort edges based on edge weight
	qsort(*sorted_edges, numedges, sizeof(float[3]), CompareArrays);

    // make sets
    int* parents = (int*) malloc(numpoints * sizeof(int));
    int* ranks = (int*) malloc(numpoints * sizeof(int));
    for (int i = 0; i < numpoints; ++i) {
        parents[i] = i;
        ranks[i] = 1;
    }

    // iterate through sorted edges and construct minimum spanning tree
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

// assigns an initialized graph random edge weights and filters out edges to help with memory/speed issues
float** generate_complete_graph(int numpoints, unsigned int* numedges) {
	float** graph = init_graph(numpoints);

	// fill graph with random numbers
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < i; ++j) {
			graph[i][j] = rand_float();

			// for numpoints > 4096 filter out edges greater than 0.03
			if (numpoints > 4096 && graph[i][j] > 0.03) {
				graph[i][j] = 0.0;
				*numedges = *numedges - 1;
			}
		}
	}

	return graph;
}

// calculates the distance between two points in a given dimensional space
float distance(float* point1, float* point2, int d) {
	double sum = 0;
	for (int i = 0; i < d; ++i) {
		sum += pow((double) point1[i] - point2[i], (double) 2);
	}
	return pow(sum, 0.5);
}

/* assigns an initialized graph with random edge weights according to dimension 
   and filters out edges we're certain won't make it into our minimum spanning tree */
float** generate_euclidean_graph(int numpoints, int dim, unsigned int* numedges) {
	float** graph = init_graph(numpoints);

	// temporary storage for coordinates
	float coordinates[numpoints][dim];
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < dim; ++j) {
			coordinates[i][j] = rand_float();
		}
	}

	// max_val of edge to keep for n > 4096 (filter out rest)
	// indexed by dim
	float max_vals[5] = {0, 0, 0.1, 0.25, 0.4};

	// fill graph with random euclidean distances
	// by randomly generating coordinates
	for (int i = 0; i < numpoints; ++i) {
		for (int j = 0; j < i; ++j) {
			graph[i][j] = distance(coordinates[i], coordinates[j], dim);

            // filters out edges based on dimensions
			if (numpoints > 4096 && graph[i][j] > max_vals[dim]) {
				graph[i][j] = 0.0;
				*numedges = *numedges - 1;
			}
		}
	}

	return graph;
}

// largely for debugging and testing purposes, prints all vertexes and edge weights in graph
void print_graph(float** graph, int numpoints) {
	// prints out adjacency matrix
    for (int i = 0; i < numpoints; ++i) {
    	for (int j = 0; j < numpoints; ++j) {
    		printf("%6.4f ", graph[i][j]);
    	}
    	printf("\n");
    }
}

/* adds up all the edge weights in a minimum spanning tree;
   also kept track of max edge weight to figure out what 
   value to filter edge weights on */
float find_mst_weight(float** graph, int numpoints, float* max_edge_weight) {
    float total_weight = 0.0;
    for (int i = 0; i < numpoints; ++i) {
        for (int j = 0; j < i; ++j) {
            total_weight += graph[i][j];

            if (graph[i][j] > *max_edge_weight) {
                *max_edge_weight = graph[i][j];
            }
        }
    }

    return total_weight;
}

int main(int argc, char *argv[]) {
    // seed rand() w/ current time
    // time_t t;
    // srand((unsigned) time(&t));

    int numpoints = atoi(argv[2]);
    int iterations = atoi(argv[3]);
    int dimension = atoi(argv[4]);

    // display entered command
    printf("%s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4]);

    float average_weight, average_time;

    /* used to keep track of the max edge weight ever found in multiple iterations of finding 
       an mst for a particular graph type. we can these use this value to filter our edges */
    float* max_edge_weight = malloc(sizeof(float));
    *max_edge_weight = 0.0;

    /* used to keep track of how many edges have been filtered out so we know how much 
       space to malloc for sorted_edges */
   	unsigned int* numedges = malloc(sizeof(unsigned int));

    float total_mst_weights = 0.0;
    float total_time = 0.0;
    time_t t;

    for (int i = 0; i < iterations; ++i) {
        srand((unsigned) time(&t));
        // need to reinitialize in every iteration because we filtered edges
        *numedges = ((unsigned int) (numpoints) * (unsigned int) (numpoints - 1))/2;

        // generate graph based on dimension
        float** graph;
        if (dimension == 0) 
            graph = generate_complete_graph(numpoints, numedges);
        else
            graph = generate_euclidean_graph(numpoints, dimension, numedges);

        // time how long it takes to find the minimum spanning tree
        clock_t start = clock(); 
        float** mst = find_mst(graph, numpoints, numedges);
        total_time += (float) (clock() - start) / CLOCKS_PER_SEC;

        // increment total_mst_weights by the weight of the mst in this iteration
        total_mst_weights += find_mst_weight(mst, numpoints, max_edge_weight);

        // free mst; the original generated graph is already freed in the find_mst function
        free_graph(mst, numpoints);
    }

    // calculate average weight and time
    average_weight = total_mst_weights / (float) iterations;
    average_time = total_time / (float) iterations;

    printf("Num edges (last iteration): %u\n", *numedges);
    printf("Num edges filtered out (last iteration): %u\n", ((numpoints) * (numpoints - 1))/2 - *numedges);
    printf("Max edge weight: %f\n", *max_edge_weight);
    printf("Average weight: %f\n", average_weight);
    printf("Average time: %f\n", average_time);

    printf("==========\n");
    return 0;
}

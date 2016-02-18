#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

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

int main(int argc, char *argv[]) {
    // seed rand() w/ current time
    time_t t;
    srand((unsigned) time(&t));

    int numpoints = atoi(argv[2]);
    int dimension = atoi(argv[4]);

    float** complete_graph = generate_complete_graph(numpoints);
    float** euc_graph = generate_euclidean_graph(numpoints, dimension);

   	// print_graph(euc_graph, numpoints);

    free_graph(complete_graph, numpoints);
    free_graph(euc_graph, numpoints);

    return 0;
}

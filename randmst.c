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

// merge and mergeSort borrowed from http://geeksquiz.com/merge-sort/ with modifications
/* Function to merge the two haves arr[l..m] and arr[m+1..r] of array arr[] */
void merge(int arr[][3], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 =  r - m;
 
    /* create temp arrays */
    int L[n1][3], R[n2][3];
 
    /* Copy data to temp arrays L[] and R[] */
    for(i = 0; i < n1; i++) {
    	L[i][0] = arr[i][0];
    	L[i][1] = arr[i][1];
    	L[i][2] = arr[i][2];
    }
    for(j = 0; j < n2; j++) {
        R[j][0] = arr[i][0];
        R[j][1] = arr[i][1];
        R[j][2] = arr[i][2];
    }
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i][0] <= R[j][0])
        {
            arr[k][0] = L[i][0];
            arr[k][1] = L[i][1];
            arr[k][2] = L[i][2];
            i++;
        }
        else
        {
            arr[k][0] = R[j][0];
            arr[k][1] = R[j][1];
            arr[k][2] = R[j][2];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there are any */
    while (i < n1)
    {
        arr[k][0] = L[i][0];
        arr[k][1] = L[i][1];
        arr[k][2] = L[i][2];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there are any */
    while (j < n2)
    {
        arr[k][0] = R[j][0];
        arr[k][1] = R[j][1];
        arr[k][2] = R[j][2];
        j++;
        k++;
    }
}
 
/* l is for left index and r is right index of the sub-array
  of arr to be sorted */
void mergeSort(int arr[][3], int l, int r)
{
    if (l < r)
    {
        int m = l+(r-l)/2; //Same as (l+r)/2, but avoids overflow for large l and h
        mergeSort(arr, l, m);
        mergeSort(arr, m+1, r);
        merge(arr, l, m, r);
   	}
}

float** find_mst(float** graph, int numpoints) {
	float** mst = init_graph(numpoints);

	int numedges = ((numpoints) * (numpoints - 1))/2;
	
	int sorted_edges[numedges][3];
	
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
	mergeSort(sorted_edges, 0, numedges - 1);

	int* parents = (int*) malloc(numpoints * sizeof(int));
	int* ranks = (int*) malloc(numpoints * sizeof(int));

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

int main(int argc, char *argv[]) {
    // // seed rand() w/ current time
    // time_t t;
    // srand((unsigned) time(&t));

    // int numpoints = atoi(argv[2]);
    // int dimension = atoi(argv[4]);

    // float** complete_graph = generate_complete_graph(numpoints);
    // float** euc_graph = generate_euclidean_graph(numpoints, dimension);

   	// // print_graph(euc_graph, numpoints);

    // free_graph(complete_graph, numpoints);
    // free_graph(euc_graph, numpoints);

    int temp[6][3];
    temp[0][0] = 3;
    temp[0][1] = 1;
    temp[0][2] = 0;
    temp[1][0] = 7;
    temp[1][1] = 2;
    temp[1][2] = 0;
    temp[2][0] = 2;
    temp[2][1] = 2;
    temp[2][2] = 1;
    temp[3][0] = 9;
    temp[3][1] = 3;
    temp[3][2] = 0;
    temp[4][0] = 1;
    temp[4][1] = 3;
    temp[4][2] = 1;
    temp[5][0] = 4;
    temp[5][1] = 3;
    temp[5][2] = 2;

    mergeSort(temp, 0, 5);

    for (int i = 0; i < 6; ++i) {
    	printf("%d %d %d\n", temp[i][0], temp[i][1], temp[i][2]);
    }

    return 0;
}

#include<stdio.h>
#include<stdlib.h>

int main(){

    int data[3] = data(3);

    for(int i = 0; i < 3; i++){

        printf("%d\n", data[i])

    }

    return 0;
    
}


//This function returns random numbers.
int data[] (int npoints){

    int a[npoints];
    
    for(int i = 0; i < npoints; i++){

        a[i] = rand();

    }

    return a;

}

int kdTree(
    int dim,
    int ndata,
    double data[],
    int kk,
    int cluster_start[],
    int cluster_size[],
    double cluster_boundary[][],
    double cluster_centroid[][],
    int cluster_assign[]
);


int bipartition();
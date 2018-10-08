#include<stdio.h>
#include<stdlib.h>
#include <time.h> 
#include <math.h>

//Number of dimesions (we define this).
#define DIM 2

//Number of points.
#define NDATA  8

//Number of Clusters to have at the end.
#define KK 4

//Do we want the same numbers.
#define SAME_NUMBERS 0

//Global counter for cluster assign.
int global_counter = 1;

//Data
int *bipartion (float data[], int i0, int im, int cluster_assign[]){

    //Get the subarray which is going to be the one to part.
    int SIZE = (im+1) - i0;
    float sub_data[SIZE];
    int sub_cluster_assign[SIZE];

    //Copy what we need to subdata.
    int counter = i0;
    for(int i = 0; i < SIZE; i++){

        sub_data[i] = data[counter];
        sub_cluster_assign[i] = cluster_assign[counter];
        counter++;
    }

    //After we get the clusters we just replace the indexes in cluster assign.

    float centroid[DIM];
    float variance[DIM];

    //Init centroid.
    for(int i = 0; i < DIM; i++){
        centroid[i] = 0;
    }
    
    //Init variance.
    for(int i = 0; i < DIM; i++){
        variance[i] = 0;
    }

    //Get centroid.
    int i = 0;
    int j = 0;

    while(i < (SIZE)){

        j = 0;

        while(j < DIM){

            centroid[j] = centroid[j] + sub_data[i];
            j++;
            i++;

        }   

    }

    //Set centroid.
    for(int i = 0; i < DIM; i++){
        centroid[i] = centroid[i] / (SIZE/DIM);
    }

    //Now we have the centroid in centroid array.

    //Get variance
    i = 0;

    while(i < (SIZE)){

        j = 0;

        while(j < DIM){

            //Find in variance (x - x_)
            variance[j] += (sub_data[i] - centroid[j]) * (sub_data[i] - centroid[j]); 
            j++;
            i++;

        }   

    }

    //Set Variance
    for(int i = 0; i < DIM; i++){
        variance[i] = variance[i] / (SIZE/DIM);
        //printf("%.2f ", centroid[i]);
    }


    //Find Max variance.
    int max_variance_index = 0;
    float current_max_variance = variance[0];

    for(int i = 0; i < DIM; i++){

        //To find max variance. (Where we are going to part the array.)
        if(current_max_variance < variance[i]){
            current_max_variance = variance[i];
            max_variance_index = i;
        }

    }
    
    //Assign the cluster to all respective numbers.
    int new_cluster_assign[SIZE];
    int tmp_cluster_assign[SIZE/DIM];

    //Initialize new_cluster_assign.
    for(int i = 0; i<(SIZE); i++){
        new_cluster_assign[i] = 0;
    }  

    //Initialize temp.
    for(int i = 0; i<(SIZE/DIM); i++){
        tmp_cluster_assign[i] = 0;
    }  
    
    //Get the highest on cluster assign. We need NDATA*DIM since we are on cluster_assign.
    int m = cluster_assign[0];
    for(int i = 0; i < NDATA*DIM; i++){
        if(m < cluster_assign[i]){
            m = cluster_assign[i];
        }
    }

    int k = 0;

    //Assign to particular cluster.
    int partition_dimension = max_variance_index;
    i = partition_dimension;

    while(i < (SIZE)){

        
        //Put a cluster based in the 
        if(sub_data[i] >= centroid[partition_dimension]){ //What if it is the same?        
                tmp_cluster_assign[k] = m+2;
        }else{
                tmp_cluster_assign[k] = m+1;
        }
        
        i += DIM;
        k++;

    }

    i = 0;
    int x = 0;

    while(i < (SIZE)){

        j = 0;

        while(j < DIM){

            new_cluster_assign[i] = tmp_cluster_assign[x];
            j++;
            i++;
        }

        x++;
    } 

    //sub_cluster_assign = new_cluster_assign;

    for(int i = 0; i<(SIZE); i++){
        sub_cluster_assign[i] = new_cluster_assign[i];
    }

    //Copy sub_cluster_assign into the indexes of cluster_assign.
    counter = i0;
    for(int i = 0; i < SIZE; i++){
        cluster_assign[counter] = sub_cluster_assign[i];
        counter++;
    }

    return cluster_assign;


}

float *search_kdtree(float data[], 
                     float query[],
                     int cluster_start[], 
                     int cluster_size[],
                     float temp_boundary[KK][2*DIM]){

    float cluster_boundaries[KK][2*DIM];
    float distances_to_clusters[KK];
    static float return_point[DIM]; //Closest point to point in cluster.

    float M;
    float m;

    int complete = 0;

    //Set cluster_boundaries
    for(int i = 0; i<KK; i++){
        for(int j = 0; j<2*DIM; j++){
            cluster_boundaries[i][j] = temp_boundary[i][j];
        }
    }

    //Find the distances to all clusters, using cluster boundaries.
    float current_sum = 0;
    for(int i = 0; i < KK; i++){

        current_sum = 0;
        for(int j=0; j<DIM; j++){

            m = cluster_boundaries[i][2*j]; //mins
            M = cluster_boundaries[i][(2*j)+1]; //max

            if(query[j] >= m && query[j] <= M){
                current_sum += 0;
            }else if(query[j] > M){ 
                current_sum += ((query[j] - M) * (query[j] - M));
            }else if(query[j] < m){
                current_sum += ((m-query[j]) * (m-query[j]));
            }else{
                printf("WTF!!!");
            }
        }

        distances_to_clusters[i] = sqrt(current_sum);

    }

    //Get the cluster with the minimun distance to point.
    int min_cluster_index = 0;
    int min_distance_cluster = distances_to_clusters[0];
    for(int i = 0; i < KK; i++){

        if(distances_to_clusters[i] < min_distance_cluster){
            min_distance_cluster = distances_to_clusters[i];
            min_cluster_index = i;
        }

    }

    printf("Distance to clusters: ");
    for(int i = 0; i < KK; i++){
        printf("%.0f ", distances_to_clusters[i]);
    }
    printf("\n");
        
    //printf("%d\n", min_cluster_index);
    
    while(complete == 0){

        //Exhaustive search through all the points in that cluster, and find DMIN
        int x = 0;
        float dmin = 0;

        float current_point[DIM];

        //Check that we are not checking the same point.
        int same = 0;

        for(int i = cluster_start[min_cluster_index]; i<cluster_start[min_cluster_index]+cluster_size[min_cluster_index]; i+=DIM){

            current_sum = 0;
            x = 0;

            for(int j = i; j<(i+DIM); j++){

                //What happens if the same point we are searching is in the cluster?
                if(query[x] == data[j]){
                    same++;
                }

                current_sum += (query[x] - data[j]) * (query[x] - data[j]);
                current_point[x] = data[j];
                printf("%.0f ", data[j]);
                x++;
            }

            //Here we know that the query point is in the cluster we are looking.
            //printf("%d\n", same);

            if(same == DIM){
                printf("Query point found on cluster %d\n", min_cluster_index);
                return query;
            }

            printf("%.2f %.2f\n", sqrt(current_sum), dmin);

            if(dmin == 0){
                dmin = sqrt(current_sum);
                for(int i = 0; i<DIM; i++){
                     return_point[i] = current_point[i];
                }
            }else{
                if(dmin > sqrt(current_sum)){
                    dmin = sqrt(current_sum);
                    for(int i = 0; i<DIM; i++){
                        return_point[i] = current_point[i];
                    }
                }
            }

        }

        //Compare dmin to all cluster distances, do an exception with min_cluster_index.
        for(int i = 0; i < KK; i++){
            if(dmin > distances_to_clusters[i] && i != min_cluster_index){
                min_cluster_index = i;
            }else{
                if(i == KK-1){
                    complete = 1;
                }
            }
        }

    }   

    return return_point;
    
}

//Also Kd-tree.
int main(void){

    //Array of numbers.
    float data[DIM*NDATA];

    //Which cluster is point at index.
    int cluster_assign[DIM*NDATA];

    //Control of how many clusters we have.
    static int nClusters = 1;
    
    //Initializing cluster_assign.
    for(int i = 0; i < (DIM*NDATA); i++){
        cluster_assign[i] = 0;
    }

    // Use current time as seed for random generator 
    if(SAME_NUMBERS == 1){
        srand(time(0));
    }
        
    
    
    //Generate random numbers.
    for(int i = 0; i < (DIM*NDATA); i++){

        data[i] =  (rand() % (100 - 0 + 1)) + 0;
        //printf("%.1f ", data[i]);
    }

    int temp_clusters[NDATA*DIM];

    //temporary store to the copy into temp_clusters
    int *tmp = (int *) malloc(DIM*NDATA);

    //Variables for sorting.
    int temp_sort;
    int c;
    int l;

    //Get each cluster size.
    int cluster_size[KK];

    //Get each cluster start. (Indexes)
    int cluster_start[KK];

    //Get cluster boundary
    float cluster_boundary[KK][2*DIM];

    int temp_counter = 0;
    int next;

    while(nClusters < KK){

        //Set cluster size and cluster start.
        for(int i = 0; i<KK; i++){
            cluster_size[i] = 0;
            cluster_start[i] = 0;
        }

        //Number of clusters that we currently.
        int variations = 1;
        c = 0;

        //Get variations (different clusters)
        while(c < (DIM*NDATA)){

                if((c + 1) != (DIM*NDATA)){
                    if(cluster_assign[c] != cluster_assign[c+1]){
                        variations++;
                    }
                }

                c++;

        }

            cluster_start[0] = 0;

            c = 0;
            l= 1;

            //Get cluster start.
            while(l < variations){

                if((c + 1) != (DIM*NDATA)){

                    if(cluster_assign[c] != cluster_assign[c+1]){
                        cluster_start[l]=c+1;
                        l++;

                    }
                }

                c++;

            }

            //Get each cluster size.
            for(int i = 0; i < variations; i++){

                //Get each size cluster. 
                for(int j = 0; j < NDATA*DIM; j+=DIM){

                    if(cluster_assign[j] == i){
                        cluster_size[i]++;
                    }

                }
            
            }

            //Size of the current first cluster after sorting.
            int first_cluster_assigns = cluster_assign[0];
            int temp_size_cluster_assign = 0;

            int i = 0;

            while(i < NDATA*DIM){

                if(first_cluster_assigns ==  cluster_assign[i]){
                    i++;
                }else{
                    temp_size_cluster_assign = i;
                    break;
                }

            }

            if(cluster_assign[0] == 0){
                temp_size_cluster_assign = NDATA*DIM;
            }
            
            /*
            *BIPARTION CALL
            */

            //Get new cluster assign in temp.
            tmp =  bipartion(data, 0, temp_size_cluster_assign-1, cluster_assign);

            //copy to temp_clusters
            for(int i = 0; i<(NDATA*DIM); i++){
                temp_clusters[i] = tmp[i];
                //printf("%d ", temp_clusters[i]);
            }

            //Sort Clusters and points.
            for(int i = 0; i < (NDATA*DIM); i+=DIM){
                for(int j = 0; j < (NDATA*DIM); j+=DIM){

                    if(temp_clusters[i] < temp_clusters[j]){
                    
                        c = j;
                        l = i;

                        for(int k = 0; k < DIM; k++){
                       
                            temp_sort = temp_clusters[c];
                            temp_clusters[c] = temp_clusters[l];
                            temp_clusters[l] = temp_sort;

                            temp_sort = data[c];
                            data[c] = data[l];
                            data[l] = temp_sort;

                            l++;
                            c++;
                        }

                    }
                }
            }
        
            //Increment by
            nClusters+=1;

            //Put into cluster_assign.
            for(int i = 0; i < DIM*NDATA; i++){
                cluster_assign[i] = temp_clusters[i];
            }

            //printf("\n");

        }//END OF MAIN WHILE: TO CALL BIPARTION.


        /*
        *Custer Size and Cluster Start
        */

        //Set cluster size and cluster start.
        for(int i = 0; i<KK; i++){
            cluster_size[i] = 0;
            cluster_start[i] = 0;
        }

        cluster_start[0] = 0;

        c = 0;
        l= 1;

        //Get cluster start.
        while(l < KK){

            if((c + 1) != (DIM*NDATA)){

                if(cluster_assign[c] != cluster_assign[c+1]){
                    cluster_start[l]=c+1;
                    l++;

                }
            }

            c++;

        }

        //Get each cluster size.
        l = 0;
        c = 0;
        int size_counter = 0;
        int temp_ = cluster_assign[0];

        while(l < KK){

            if(cluster_assign[c] == temp_){

                size_counter++;
                c++;

            }else{
                cluster_size[l] = size_counter;
                l++;
                size_counter = 0;
                temp_ = cluster_assign[c];
            }

            

        }

        //Print Cluster Assign, Data
        printf("--------------\n");

        printf("Cluster Assign: ");
        for(int i = 0; i < DIM*NDATA; i++){
            printf("%d ", cluster_assign[i]);
        }

        printf("\n--------------\n");

        printf("Data: ");
        for(int i = 0; i < DIM*NDATA; i++){
            printf("%.0f ", data[i]);
        }

        printf("\n--------------\n");

        printf("Cluster Size: ");
        for(int i = 0; i < KK; i++){
            cluster_size[i] = cluster_size[i];
            printf("%d ", cluster_size[i]);
        }

        printf("\n--------------\n");

        printf("Cluster Start: ");
        for(int i = 0; i < KK; i++){
            printf("%d ", cluster_start[i]);
        }

        printf("\n--------------\n");

        /*
            To get cluster boundary.
        */

        int counter= 0;
        int k ;
        int L;
        int size_c;
        
        while(counter < KK){   

            k=0;
            L= cluster_start[counter];
            size_c = cluster_size[counter];

            for(int j =L;j<(L+size_c);){

                if(j==L){

                    while(k<2*DIM){

                        cluster_boundary[counter][k]=data[j];
                        k++;

                        cluster_boundary[counter][k]=data[j];
                        k++;
                        j++;

                    }

                }

            else{

                k=0;
                while(k<2*DIM){

                    if(cluster_boundary[counter][k]>data[j])
                        cluster_boundary[counter][k]=data[j];   
                    k++;
                    if(cluster_boundary[counter][k]<data[j])
                        cluster_boundary[counter][k]=data[j];
                    k++;

                    j++;

                    }

                }
            }
            counter++;
        }

        printf("Cluster Boundaries: \n"); 
        //Print cluster boundaries.
        for(int i = 0; i < KK; i++){

            for(int j = 0; j<2*DIM; j++){

                printf("%f ", cluster_boundary[i][j]);

            }

            printf("\n");

        }
        printf("-----------------\n");

        //Call search

        float query[DIM];
        query[0] = 65;
        query[1] = 52;

        printf("Query Point: ");
        for(int i = 0; i<DIM; i++){
            printf("%.0f ", query[i]);
        }
        printf("\n");

        float *rpoint;

        rpoint = search_kdtree(data, query, cluster_start, cluster_size, cluster_boundary);
        //printf("here2");

        float end_point[DIM];

        //Print return point.
       printf("Return Point: ");
        for(int i = 0; i<DIM; i++){
            end_point[i] = rpoint[i];
            printf("%.0f ", end_point[i]);
        }
        printf("\n--------------\n");
        
        return 0;

}



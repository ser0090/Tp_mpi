#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 8
#define MASTER 0
#define TAG_A 1
#define TAG_B 2
#define TAG_C 3

int main(int argc, char** argv) {

    int world_size;
    int rank;

    int matrix_a[SIZE][SIZE];
    int matrix_b[SIZE][SIZE];
    int matrix_c[SIZE][SIZE];

    int sub_a[SIZE/2][SIZE];
    int sub_b[SIZE][SIZE/2];
    int sub_c[SIZE/2][SIZE/2];

    MPI_Init(NULL, NULL);   // Initialize the MPI environment
    MPI_Status stat;        // required variable for receive routines

    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // Get the rank of the process

    if(rank == MASTER){
        for (int i = 0; i < SIZE; ++i) {
            for (int j = 0; j < SIZE; ++j) {
                matrix_a[i][j]=j+i*SIZE;
                matrix_b[i][j]=0;
                //printf("%d ",matrix_a[i][j]);
            }
            //printf("\n");
        }
        for (int i = 0; i < SIZE; ++i)
            matrix_b[i][i]=1; //matriz identidad

    }

    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < 2 ; ++i) {
            for (int j = 0; j <SIZE/2 ; j++) {
                if (rank == MASTER){
                    MPI_Send(&matrix_a[j+k*SIZE/2][i*SIZE/2], SIZE/2, MPI_INT, (i+2*k), TAG_A, MPI_COMM_WORLD);
                    MPI_Send(&matrix_b[j+k*SIZE/2][i*SIZE/2], SIZE/2, MPI_INT, (i+2*k), TAG_B, MPI_COMM_WORLD);
                }
                if(rank == (i+2*k)){
                    MPI_Recv(&sub_a[j][(rank&1)*SIZE/2], SIZE/2, MPI_INT, MASTER, TAG_A, MPI_COMM_WORLD, &stat);
                    //TODO: revisar j-(rank&1) esta medio hardcode
                    MPI_Recv(&sub_b[j-(rank&1)+k*SIZE/2][i*SIZE/2], SIZE/2, MPI_INT, MASTER, TAG_B, MPI_COMM_WORLD, &stat);
                }
            }
        }
    }


    if(!(rank&1)){ //is rank es par para A
        for (int j = 0; j < SIZE/2; ++j) {
            MPI_Send(sub_a[j],SIZE/2, MPI_INT, rank+1, TAG_A, MPI_COMM_WORLD);
            MPI_Recv(&sub_a[j][SIZE/2], SIZE/2, MPI_INT, rank+1, TAG_A, MPI_COMM_WORLD, &stat);
        }
    }
    else{
        for (int j = 0; j < SIZE/2; ++j) {
            MPI_Recv(sub_a[j], SIZE/2, MPI_INT, rank-1, TAG_A, MPI_COMM_WORLD, &stat);
            MPI_Send(&sub_a[j][SIZE/2], SIZE/2, MPI_INT, rank-1, TAG_A, MPI_COMM_WORLD);
        }
    }

    if(!(rank/2)){ // Pasaje de Matrix B
        for (int j = 0; j < SIZE/2; ++j) {
            MPI_Send(sub_b[j],SIZE/2, MPI_INT, rank+2, TAG_B, MPI_COMM_WORLD);
            MPI_Recv(sub_b[j+SIZE/2], SIZE/2, MPI_INT, rank+2, TAG_B, MPI_COMM_WORLD, &stat);
        }
    }
    else{
        for (int j = 0; j < SIZE/2; ++j) {
            MPI_Recv(sub_b[j], SIZE/2, MPI_INT, rank-2, TAG_B, MPI_COMM_WORLD, &stat);
            MPI_Send(sub_b[j+SIZE/2], SIZE/2, MPI_INT, rank-2, TAG_B, MPI_COMM_WORLD);
        }
    }

    int tmp;
    // PROUCTO dejando fija la columna de B
    for (int i=0; i<SIZE/2; i++){ //i para las filas de la matriz resultante
        for (int j=0; j<SIZE/2 ;j++){ // i para las columnas de la matriz resultante
            tmp = 0 ;
            for (int k=0; k<SIZE; k++){ //k para realizar la multiplicacion de los elementos
                tmp += sub_a[j][k] * sub_b[k][i];
            }
            sub_c[j][i] = tmp;
        }
    }

    if(rank == MASTER) {
        for (int j=0; j<SIZE/2; ++j) {
            MPI_Recv(&matrix_c[j][SIZE/2], SIZE/2, MPI_INT, 1, TAG_C, MPI_COMM_WORLD, &stat);
            MPI_Recv(matrix_c[j+SIZE/2], SIZE/2, MPI_INT, 2, TAG_C, MPI_COMM_WORLD, &stat);
            MPI_Recv(&matrix_c[j+SIZE/2][SIZE/2], SIZE/2, MPI_INT, 3, TAG_C, MPI_COMM_WORLD, &stat);

            MPI_Send(sub_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD);
            MPI_Recv(matrix_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD, &stat);
        }
        //IMPRIME C
        for (int i = 0; i < SIZE; ++i) {
            for (int j = 0; j < SIZE; ++j) {
                printf("%d ", matrix_c[i][j]);
            }
            printf("\n");
        }
    }
    else {
        for (int j = 0; j < SIZE/2; ++j) {
            MPI_Send(sub_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

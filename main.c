#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define SIZE 2000
#define HALF_SIZE SIZE/2
#define MASTER 0
#define TAG_A 1
#define TAG_B 2
#define TAG_C 3
#define ERROR 1
#define SUCCESS 0

#define COLOR_GREEN   "\x1b[1;32m"
#define COLOR_RESET   "\x1b[0m"


typedef int32_t my_matrix_t [SIZE][SIZE];
typedef int32_t my_sub_a_t [HALF_SIZE][SIZE];
typedef int32_t my_sub_b_t [SIZE][HALF_SIZE];
typedef int32_t my_sub_c_t [HALF_SIZE][HALF_SIZE];

int control(int32_t **matrix, int32_t size);

int main(int argc, char** argv) {

    int32_t rank;
    register int32_t tmp;

    double start_time;
    double start;
    my_matrix_t* matrix_a;
    my_matrix_t* matrix_b;
    my_matrix_t* matrix_c;

    my_sub_c_t* sub_c;
    my_sub_a_t* sub_a;
    my_sub_b_t* sub_b;

    MPI_Init(NULL, NULL);   // Initialize the MPI environment

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // Get the rank of the process

    /**#################################################################
    ######################  INICIALIZACION  ###########################
    ##################################################################**/

    sub_a = (my_sub_a_t*) malloc((HALF_SIZE)*SIZE*sizeof(int32_t)+4);
    sub_b = (my_sub_b_t*) malloc(SIZE*(HALF_SIZE)*sizeof(int32_t));
    sub_c = (my_sub_c_t*) malloc((HALF_SIZE)*(HALF_SIZE)*sizeof(int32_t));

    if(rank == MASTER){

        start = omp_get_wtime();
	    start_time = start;

        printf("RANK %d:  size: %d\n", rank, SIZE);

        matrix_a = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));
        matrix_b = (my_matrix_t*) calloc(SIZE*SIZE, sizeof(int32_t));
        matrix_c = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));

        for (int32_t j = 0; j < SIZE; ++j) {
            for (int32_t i = 0; i < SIZE; ++i) {
                (*matrix_a)[i][j]=j+i*SIZE;
            }
        }
        for (int32_t i = 0; i < SIZE; ++i) {
            (*matrix_b)[i][i] = 1; //matriz identidad
        }
        // el master copia los valores a la sub_matrix
        for (int32_t j = 0; j <HALF_SIZE ; j++) {
            for (int32_t i = 0; i < HALF_SIZE ; ++i) {
                (*sub_a)[i][j] = (*matrix_a)[i][j];
                (*sub_b)[i][j] = (*matrix_b)[i][j];
            }
        }

        printf("TIME: iniciar la matriz: %f\n",omp_get_wtime()-start);

        /**#################################################################
        #################  FRACCIONAMIENTO DE MATRICES  ###################
        ##################################################################**/

        start = omp_get_wtime();
        for (int32_t j = 0; j < HALF_SIZE; j++) {
            MPI_Send(&((*matrix_a)[j][HALF_SIZE]), HALF_SIZE, MPI_INT32_T,1, TAG_A, MPI_COMM_WORLD);
            MPI_Send(&((*matrix_b)[j][HALF_SIZE]), HALF_SIZE, MPI_INT32_T,1, TAG_B, MPI_COMM_WORLD);

            MPI_Send(&((*matrix_a)[j + HALF_SIZE][0]), HALF_SIZE, MPI_INT32_T, 2, TAG_A, MPI_COMM_WORLD);
            MPI_Send(&((*matrix_b)[j + HALF_SIZE][0]), HALF_SIZE, MPI_INT32_T, 2, TAG_B, MPI_COMM_WORLD);

            MPI_Send(&((*matrix_a)[j + HALF_SIZE][HALF_SIZE]), HALF_SIZE, MPI_INT32_T, 3, TAG_A, MPI_COMM_WORLD);
            MPI_Send(&((*matrix_b)[j + HALF_SIZE][HALF_SIZE]), HALF_SIZE, MPI_INT32_T, 3, TAG_B, MPI_COMM_WORLD);
        }
        printf("TIME: fraccionar: %f\n",omp_get_wtime()-start);
    }
    else{
        for (int32_t j = 0; j <HALF_SIZE ; j++) {
            MPI_Recv(&((*sub_a)[j][(rank & 1) * HALF_SIZE]), HALF_SIZE, MPI_INT32_T, MASTER, TAG_A, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
            MPI_Recv((*sub_b)[j + (rank > 2) * HALF_SIZE], HALF_SIZE, MPI_INT32_T, MASTER, TAG_B, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
        }
    }

    /**#################################################################
    #######################  ENVIO DE SUB MATRICES  ###################
    ##################################################################**/

    if(!(rank & 1)){ //is rank es par para A
        for (int32_t j = 0; j < HALF_SIZE; ++j) {
            MPI_Sendrecv((*sub_a)+j, HALF_SIZE, MPI_INT32_T, rank+1, TAG_A,&((*sub_a)[j][HALF_SIZE]),HALF_SIZE, MPI_INT32_T,
                    rank+1, TAG_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else{
        for (int32_t j = 0; j < HALF_SIZE; ++j) {
            MPI_Sendrecv(&((*sub_a)[j][HALF_SIZE]), HALF_SIZE, MPI_INT32_T, rank-1, TAG_A, (*sub_a)+j, HALF_SIZE, MPI_INT32_T,
                    rank-1,TAG_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if(!(rank >> 1)){ // Pasaje de Matrix B  (rank >> 1) -> rank/2
        for (int32_t j = 0; j < HALF_SIZE; ++j) {
            MPI_Sendrecv((*sub_b)+j, HALF_SIZE, MPI_INT32_T, rank+2, TAG_B, (*sub_b)+j+HALF_SIZE, HALF_SIZE, MPI_INT32_T,
                    rank+2, TAG_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else{
        for (int32_t j = 0; j < HALF_SIZE; ++j) {
            MPI_Sendrecv((*sub_b)+j+HALF_SIZE, HALF_SIZE, MPI_INT32_T, rank-2, TAG_B, (*sub_b)+j, HALF_SIZE, MPI_INT32_T,
                    rank-2,TAG_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    printf("RANK %d: fin de subenvios\n", rank);

    /**#################################################################
    #######################   PROCESAMIENTO   #########################
    ##################################################################**/

    // PROUCTO dejando fija la columna de B
    for (int32_t i=0; i<HALF_SIZE; i++){ //i para las filas de la matriz resultante
        for (int32_t j=0; j<HALF_SIZE ;j++){ // i para las columnas de la matriz resultante
            tmp = 0 ;
            for (int32_t k=0; k<SIZE; k++){ //k para realizar la multiplicacion de los elementos
                tmp += (*sub_a)[j][k] * (*sub_b)[k][i];
            }
            (*sub_c)[j][i] = tmp;
        }
    }
    printf("RANK %d: procesar todo\n",rank);

    /**#################################################################
    ############################ RESULTADO   #########################
    ##################################################################**/

    if(rank == MASTER) {
        start = omp_get_wtime();
        for (int32_t j=0; j<HALF_SIZE; ++j) {
            MPI_Recv(&((*matrix_c)[j][HALF_SIZE]), HALF_SIZE, MPI_INT32_T, 1, TAG_C, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv((*matrix_c)+j+HALF_SIZE, HALF_SIZE, MPI_INT32_T, 2, TAG_C, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&((*matrix_c)[j+HALF_SIZE][HALF_SIZE]), HALF_SIZE, MPI_INT32_T, 3, TAG_C, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
        }
        printf("RANK %d: Termino todo\n",rank);

        for (int32_t i = 0; i < HALF_SIZE ; ++i) {
            for (int32_t j = 0; j <HALF_SIZE ; j++) {
                (*matrix_c)[i][j] = (*sub_c)[i][j];
            }
        }
        printf("TIME: construir C: %f\n", omp_get_wtime() - start);
        printf("%sTIME: Total: %f%s\n", COLOR_GREEN, omp_get_wtime() - start_time, COLOR_RESET);
        free(matrix_a);
        free(matrix_b);
        // COMPROBACION
        if(control((int32_t**)matrix_c,SIZE*SIZE)==SUCCESS)
            printf("CHECK: Todo OK wacho!\n");
        else
            printf("CHECK: SE PUDRIO TODO!\n");

        free(matrix_c);
    }
    else {
        for (int32_t j = 0; j < HALF_SIZE; ++j) {
            MPI_Send((*sub_c)+j, HALF_SIZE, MPI_INT32_T, MASTER, TAG_C, MPI_COMM_WORLD);
        }
        printf("RANK %d: Termino todo\n",rank);
    }
    free(sub_a);
    free(sub_b);
    free(sub_c);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

int control(int32_t **matrix, int32_t size) {
    int32_t* aux= (int32_t*)matrix;
    for(int32_t i=0; i<size-1 ;i++){
        if(*(aux+i) >= *(aux+i+1)) return ERROR;
    }
    return SUCCESS;
}

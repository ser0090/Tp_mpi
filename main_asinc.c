#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define SIZE 1500
#define MASTER 0
#define TAG_A 1
#define TAG_B 2
#define TAG_C 3
#define ERROR 1
#define SUCCESS 0
#define SUB_SIZE SIZE/2

#define COLOR_GREEN   "\x1b[1;32m"
#define COLOR_RESET   "\x1b[0m"

typedef int32_t my_matrix_t [SIZE][SIZE];

typedef int32_t my_sub_a_t [SUB_SIZE][SIZE];
typedef int32_t my_sub_b_t [SIZE][SUB_SIZE];
typedef int32_t my_sub_c_t [SUB_SIZE][SUB_SIZE];

int control(int32_t **matrix, int32_t size);


int main(int argc, char** argv) {

    int32_t world_size;
    int32_t rank;
    int32_t shift;
    register int32_t tmp;

    double start_time;
    double start;

    int32_t *ptr_a;

    my_matrix_t* matrix_a;
    my_matrix_t* matrix_b;
    my_matrix_t* matrix_c;

    my_sub_c_t* sub_c;
    my_sub_a_t* sub_a;
    my_sub_b_t* sub_b;

    MPI_Request *fragment_a_req[4] = {NULL,NULL,NULL,NULL};
    MPI_Request *fragment_b_req[4] = {NULL,NULL,NULL,NULL};
    MPI_Request *fragment_c_req[4] = {NULL,NULL,NULL,NULL};

    MPI_Request *ptr_rq_a;
    MPI_Request *ptr_rq_b;
    MPI_Request *ptr_rq_c;

    MPI_Init(NULL, NULL);   // Initialize the MPI environment

    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // Get the rank of the process

    /**#################################################################
    ######################  INICIALIZACION  ###########################
    ##################################################################**/

    sub_a = (my_sub_a_t*) malloc((SUB_SIZE)*SIZE*sizeof(int32_t)+4);
    sub_b = (my_sub_b_t*) malloc(SIZE*(SUB_SIZE)*sizeof(int32_t));
    sub_c = (my_sub_c_t*) malloc((SUB_SIZE)*(SUB_SIZE)*sizeof(int32_t));

    if(rank == MASTER){

        start = omp_get_wtime();
        MPI_Request aux_req;

        start_time = start;

	    matrix_a = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));
        matrix_b = (my_matrix_t*) calloc(SIZE*SIZE, sizeof(int32_t));

        printf("RANK %d:  size: %d\n", rank, SIZE);

        for (int32_t i = 0; i < SIZE; ++i) {
            (*matrix_b)[i][i] = 1; //matriz identidad
        }

	    printf("TIME: iniciar la matriz: %f\n",omp_get_wtime()-start);

        /**#################################################################
        #################  FRACCIONAMIENTO DE MATRICES  ###################
        ##################################################################**/

        start = omp_get_wtime();

        //Envio parcial de B
        for (int32_t j = 0; j < SUB_SIZE; j++) {
            MPI_Isend(&((*matrix_b)[j+SUB_SIZE][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 1, TAG_B, MPI_COMM_WORLD, &aux_req);
            MPI_Isend(&((*matrix_b)[j][0]), SUB_SIZE, MPI_INT32_T, 2, TAG_B, MPI_COMM_WORLD, &aux_req);
            MPI_Isend(&((*matrix_b)[j+SUB_SIZE][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 3, TAG_B, MPI_COMM_WORLD, &aux_req);
            //MPI_Request_free(&aux_req);
        }

        for (int32_t i = 0; i < SIZE; ++i) {
            for (int32_t j = 0; j < SIZE; ++j) {
                (*matrix_a)[i][j] = j + i * SIZE;
            }
        }

        //envio de A
        for (int32_t j = 0; j < SUB_SIZE; j++) {
            MPI_Isend(&((*matrix_a)[j][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 1, TAG_A, MPI_COMM_WORLD, &aux_req);
            MPI_Isend(&((*matrix_a)[j + SUB_SIZE][0]), SUB_SIZE, MPI_INT32_T, 2, TAG_A, MPI_COMM_WORLD,&aux_req);
            MPI_Isend(&((*matrix_a)[j + SUB_SIZE][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 3, TAG_A, MPI_COMM_WORLD,&aux_req);
            //MPI_Request_free(&aux_req);
        }

        //envio final de B
        for (int32_t j = 0; j < SUB_SIZE; j++) {
            MPI_Isend(&((*matrix_b)[j][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 1, TAG_B, MPI_COMM_WORLD, &aux_req);
            MPI_Isend(&((*matrix_b)[j + SUB_SIZE][0]), SUB_SIZE, MPI_INT32_T, 2, TAG_B, MPI_COMM_WORLD, &aux_req);
            MPI_Isend(&((*matrix_b)[j][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 3, TAG_B, MPI_COMM_WORLD, &aux_req);
            //MPI_Request_free(&aux_req);
        }

        MPI_Request_free(&aux_req);

        // el master copia los valores a la sub_matrix
        for (int32_t i = 0; i < SUB_SIZE ; ++i) {
            for (int32_t j = 0; j <SIZE ; j++) {
                (*sub_a)[i][j] = (*matrix_a)[i][j];
            }
        }
        for (int32_t i = 0; i < SIZE ; ++i) {
            for (int32_t j = 0; j <SUB_SIZE ; j++) {
                (*sub_b)[i][j] = (*matrix_b)[i][j];
            }
        }
        printf("TIME: fraccionar: %f\n",omp_get_wtime()-start);
    }
    else{
        fragment_a_req[rank] = (MPI_Request*) malloc(SUB_SIZE * sizeof(MPI_Request));
        fragment_b_req[rank] = (MPI_Request*) malloc(SUB_SIZE * sizeof(MPI_Request));

        if(fragment_a_req[rank] == NULL || fragment_b_req[rank] == NULL){
            printf("CALLOC_ERROR\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        shift = (rank != 2)*SUB_SIZE;
        ptr_rq_a = &(fragment_b_req[rank][0]);
        for (int32_t j = 0; j <SUB_SIZE ; j++) {
            MPI_Irecv((*sub_b)[j+shift], SUB_SIZE, MPI_INT32_T, MASTER, TAG_B, MPI_COMM_WORLD, ptr_rq_a++);
        }
        shift = (rank & 1)*SUB_SIZE;
        ptr_rq_b = &(fragment_a_req[rank][0]);
        for (int32_t j = 0; j <SUB_SIZE ; j++) {
            MPI_Irecv(&((*sub_a)[j][shift]), SUB_SIZE, MPI_INT32_T, MASTER, TAG_A, MPI_COMM_WORLD, ptr_rq_b++);
        }

        MPI_Waitall(SUB_SIZE, fragment_b_req[rank], MPI_STATUSES_IGNORE);
    }

    /**#################################################################
    #######################   SEMI PROCESAMIENTO   ####################
    ##################################################################**/

    printf("RANK %d: Semiprocesamiento \n",rank);
    if (rank != MASTER) { // Pasaje de Matrix B
        shift = (!(rank & 1)) * SUB_SIZE;
        ptr_rq_b = &(fragment_b_req[rank][0]);
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            MPI_Irecv(((*sub_b)[j + shift]), SUB_SIZE, MPI_INT32_T, MASTER, TAG_B, MPI_COMM_WORLD, ptr_rq_b++);
        }
    }

    // PROUCTO dejando fija la columna de B
    //despliegue del procesamiento para esperar lo datos de a
    shift = (rank & 1) * SUB_SIZE;
    for (int32_t j = 0; j < SUB_SIZE; j++) { // i para las columnas de la matriz resultante
        ptr_a = &((*sub_a)[j][shift]);
        tmp = 0;
        if (rank != MASTER) MPI_Wait(&(fragment_a_req[rank][j]), MPI_STATUS_IGNORE);
        for (int32_t k = 0; k < SUB_SIZE; k++) { //k para realizar la multiplicacion de los elementos
            tmp += *(ptr_a++) * (*sub_b)[k + shift][0];
        }
        (*sub_c)[j][0] = tmp;
    }

    ///#################################################################
    //printf("RANK %d: Intercambio de A \n",rank);

    if (!(rank & 1)) { //is rank es par para A
        fragment_a_req[rank + 1] = (MPI_Request *) malloc(SUB_SIZE * sizeof(MPI_Request));
        MPI_Request aux_req;
        ptr_rq_a = &(fragment_a_req[rank + 1][0]);
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            //if(rank != MASTER) MPI_Wait(&(fragment_a_req[rank][j]), MPI_STATUS_IGNORE);
            MPI_Isend((*sub_a)+j, SUB_SIZE, MPI_INT32_T, rank+1, TAG_A, MPI_COMM_WORLD, &aux_req);
            //MPI_Request_free(&aux_req);
            MPI_Irecv(&((*sub_a)[j][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, rank+1, TAG_A, MPI_COMM_WORLD, ptr_rq_a++);
        }
        MPI_Request_free(&aux_req);
    } else {
        fragment_a_req[rank-1] = (MPI_Request *) malloc(SUB_SIZE * sizeof(MPI_Request));
        MPI_Request aux_req;
        ptr_rq_a = &(fragment_a_req[rank-1][0]);
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            MPI_Irecv((*sub_a)+j, SUB_SIZE, MPI_INT32_T, rank-1, TAG_A, MPI_COMM_WORLD, ptr_rq_a++);
            //if (rank != MASTER) MPI_Wait(&(fragment_a_req[rank][j]), MPI_STATUS_IGNORE);
            MPI_Isend(&((*sub_a)[j][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, rank - 1, TAG_A, MPI_COMM_WORLD, &aux_req);
            //MPI_Request_free(&aux_req);
        }
        MPI_Request_free(&aux_req);
    }

    ///##################################################################
    //printf("RANK %d: Semiprocesamiento  \n",rank);

    for (int32_t i = 1; i < SUB_SIZE; i++) { //i para las filas de la matriz resultante
        for (int32_t j = 0; j < SUB_SIZE; j++) { // i para las columnas de la matriz resultante
            tmp = 0;
            ptr_a = &((*sub_a)[j][shift]);
            for (int32_t k = 0; k < SUB_SIZE; k++) { //k para realizar la multiplicacion de los elementos
                tmp += *(ptr_a++) * (*sub_b)[k + shift][i];
            }
            (*sub_c)[j][i] = tmp;
        }
    }


    /**#################################################################
    #######################   PROCESAMIENTO   #########################
    ##################################################################**/

    // PROUCTO dejando fija la columna de B

    if(rank!=MASTER)
        MPI_Waitall(SUB_SIZE, fragment_b_req[rank], MPI_STATUSES_IGNORE);

    fragment_c_req[rank] = (MPI_Request*) malloc(SUB_SIZE * sizeof(MPI_Request));
    shift = (!(rank&1)) * SUB_SIZE;

        //despliegue del procesamiento para esperar lo datos de a
    for (int32_t j=0; j<SUB_SIZE ;j++){ // i para las columnas de la matriz resultante
        tmp = 0 ;
        ptr_a = &((*sub_a)[j][shift]);
        if (rank != MASTER) {
            MPI_Wait(&(fragment_a_req[rank][j]), MPI_STATUS_IGNORE);
        }
        for (int32_t k=0; k<SUB_SIZE; k++){ //k para realizar la multiplicacion de los elementos
            tmp += *(ptr_a++) * (*sub_b)[k+shift][0];
        }
        (*sub_c)[j][0] += tmp;
    }

    if(rank != MASTER)
        MPI_Isend((*sub_c)[0], SUB_SIZE, MPI_INT32_T, MASTER, rank, MPI_COMM_WORLD, &(fragment_c_req[rank][0]));

    for (int32_t i=1; i<SUB_SIZE; i++){ //i para las filas de la matriz resultante
        for (int32_t j=0; j<SUB_SIZE ;j++){ // i para las columnas de la matriz resultante
            tmp = 0 ;
            ptr_a = &((*sub_a)[j][shift]);
            for (int32_t k=0; k<SUB_SIZE; k++){ //k para realizar la multiplicacion de los elementos
                tmp += *(ptr_a++) * (*sub_b)[k+shift][i];
            }
            (*sub_c)[j][i] += tmp;
        }
        if(rank!=MASTER)
            MPI_Isend((*sub_c)[i], SUB_SIZE, MPI_INT32_T, MASTER, rank, MPI_COMM_WORLD, &(fragment_c_req[rank][i]));
    }


    /**#################################################################
    #######################   RECONSTRUCCION   #########################
    ##################################################################**/
    if(rank == MASTER) {
        for(int i=1;i<4;i++){
            fragment_c_req[i] = (MPI_Request*) malloc(SUB_SIZE * sizeof(MPI_Request));
        }
        start = omp_get_wtime();
        matrix_c = (my_matrix_t*) malloc(SIZE*SIZE * sizeof(int32_t));

        for (int32_t j=0; j<SUB_SIZE; ++j) {
            MPI_Irecv(&((*matrix_c)[j][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 1, 1, MPI_COMM_WORLD,
                    &(fragment_c_req[1][j]));
            MPI_Irecv(&((*matrix_c)[j+SUB_SIZE][0]), SUB_SIZE, MPI_INT32_T, 2, 2, MPI_COMM_WORLD,
                    &(fragment_c_req[2][j]));
            MPI_Irecv(&((*matrix_c)[j+SUB_SIZE][SUB_SIZE]), SUB_SIZE, MPI_INT32_T, 3, 3, MPI_COMM_WORLD,
                    &(fragment_c_req[3][j]));
        }

        for (int32_t i = 0; i < SUB_SIZE ; ++i) {
            for (int32_t j = 0; j <SUB_SIZE ; j++) {
              (*matrix_c)[i][j]=(*sub_c)[i][j];
            }
        }
        printf("RANK %d: Termino todo\n",rank);
        MPI_Waitall(SUB_SIZE, fragment_c_req[1], MPI_STATUS_IGNORE);
        MPI_Waitall(SUB_SIZE, fragment_c_req[2], MPI_STATUS_IGNORE);
        MPI_Waitall(SUB_SIZE, fragment_c_req[3], MPI_STATUS_IGNORE);

        printf("TIME: construir C: %f\n", omp_get_wtime() - start);
        printf("%sTIME: Total: %f%s\n", COLOR_GREEN, omp_get_wtime() - start_time, COLOR_RESET);

        free(matrix_a);
        free(matrix_b);

        // COMPROBACION
        if(control((int32_t**)matrix_c,SIZE*SIZE) == SUCCESS)
            printf("CHECK: Todo OK wacho ^^\n");
        else
            printf("CHECK: SE PUDRIO TODO T.T\n");

        free(matrix_c);
    }
    else {
        printf("RANK %d: Termino todo\n", rank);

        MPI_Waitall(SUB_SIZE,fragment_c_req[rank],MPI_STATUSES_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    free(sub_a);
    free(sub_b);
    free(sub_c);
    MPI_Finalize();

    exit(EXIT_SUCCESS);
}

int control(int32_t **matrix, int32_t size)
{
  int32_t* aux= (int32_t*)matrix;
  for(int32_t i=0;i<size-1;i++){
    if(*(aux+i)>=*(aux+i+1)) return ERROR;
  }
    return SUCCESS;
}

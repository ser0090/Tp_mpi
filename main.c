#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 4096
#define MASTER 0
#define TAG_A 1
#define TAG_B 2
#define TAG_C 3
#define ERROR 1
#define SUCCESS 0

typedef int32_t my_matrix_t [SIZE][SIZE];

typedef int32_t my_sub_a_t [SIZE/2][SIZE];
typedef int32_t my_sub_b_t [SIZE][SIZE/2];
typedef int32_t my_sub_c_t [SIZE/2][SIZE/2];

int control(int32_t **matrix, int32_t size);

int main(int argc, char** argv) {

    int32_t world_size;
    int32_t rank;

    my_matrix_t* matrix_a;
    my_matrix_t* matrix_b;
    my_matrix_t* matrix_c;

    my_sub_c_t* sub_c;
    my_sub_a_t* sub_a;
    my_sub_b_t* sub_b;

    //int32_t sub_a[SIZE/2][SIZE];
    //int32_t sub_c[SIZE/2][SIZE/2];
    //int32_t sub_b[SIZE][SIZE/2];

    if(rank==MASTER){
      void* aux;
      //aux =
      matrix_a = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));
      matrix_b = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));
      matrix_c = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));

    }

    sub_a = (my_sub_a_t*) malloc((SIZE/2)*SIZE*sizeof(int32_t));
    sub_b = (my_sub_b_t*) malloc(SIZE*(SIZE/2)*sizeof(int32_t));
    sub_c = (my_sub_c_t*) malloc((SIZE/2)*(SIZE/2)*sizeof(int32_t));

    //int32_t sub_a[SIZE/2][SIZE];
    //int32_t sub_c[SIZE/2][SIZE/2];
    //int32_t sub_b[SIZE][SIZE/2];
//    printf("antes del mpi init\n");
    MPI_Init(NULL, NULL);   // Initialize the MPI environment
    MPI_Status stat;        // required variable for receive routines

    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // Get the rank of the process
    printf("por hacer algo rankl %d\n", rank);
    if(rank == MASTER){
        for (int32_t i = 0; i < SIZE; ++i) {
            for (int32_t j = 0; j < SIZE; ++j) {
                (*matrix_a)[i][j]=j+i*SIZE;
                (*matrix_b)[i][j]=0;
                //printf("%d ",matrix_a[i][j]);
            }
//            printf("columna %d",i);
        }
        for (int32_t i = 0; i < SIZE; ++i)
            (*matrix_b)[i][i]=1; //matriz identidad

    }
//    printf("Estoy por enviar/recibir primer transmision rank:%d\n",rank );

    for (int32_t k = 0; k < 2; ++k) {
        for (int32_t i = 0; i < 2 ; ++i) {
            for (int32_t j = 0; j <SIZE/2 ; j++) {
                if (rank == MASTER && (i+2*k)){
//                    printf("por enviar a rank %d\n", (i+2*k));
                    MPI_Send(&((*matrix_a)[j+k*SIZE/2][i*SIZE/2]), SIZE/2, MPI_INT, (i+2*k), TAG_A, MPI_COMM_WORLD);
                    MPI_Send(&((*matrix_b)[j+k*SIZE/2][i*SIZE/2]), SIZE/2, MPI_INT, (i+2*k), TAG_B, MPI_COMM_WORLD);
//                    printf("enviado a rank %d\n", (i+2*k));
                }
                if(rank == (i+2*k) && (i+2*k)){
                    MPI_Recv(&((*sub_a)[j][(rank&1)*SIZE/2]), SIZE/2, MPI_INT, MASTER, TAG_A, MPI_COMM_WORLD, &stat);
                    //TODO: revisar j-(rank&1) esta medio hardcode
                    MPI_Recv(&((*sub_b)[j-(rank&1)+k*SIZE/2][i*SIZE/2]), SIZE/2, MPI_INT, MASTER, TAG_B, MPI_COMM_WORLD, &stat);
                }
            }
        }
    }

    if(rank==0){ // el master copia los valores a la sub_matrix
//      printf("Estoy por copiar mi sub rank:%d\n",rank );
      for (int32_t i = 0; i < SIZE/2 ; ++i) {
          for (int32_t j = 0; j <SIZE/2 ; j++) {
            (*sub_a)[i][j]=(*matrix_a)[i][j];
            (*sub_b)[i][j]=(*matrix_b)[i][j];
          }
      }
    }

//    printf("Estoy por enviar sub rank:%d\n",rank );

    if(!(rank&1)){ //is rank es par para A
        for (int32_t j = 0; j < SIZE/2; ++j) {
            MPI_Send(&((*sub_a)[j]),SIZE/2, MPI_INT, rank+1, TAG_A, MPI_COMM_WORLD);
            MPI_Recv(&((*sub_a)[j][SIZE/2]), SIZE/2, MPI_INT, rank+1, TAG_A, MPI_COMM_WORLD, &stat);
        }
    }
    else{
        for (int32_t j = 0; j < SIZE/2; ++j) {
            MPI_Recv(&((*sub_a)[j]), SIZE/2, MPI_INT, rank-1, TAG_A, MPI_COMM_WORLD, &stat);
            MPI_Send(&((*sub_a)[j][SIZE/2]), SIZE/2, MPI_INT, rank-1, TAG_A, MPI_COMM_WORLD);
        }
    }

    if(!(rank/2)){ // Pasaje de Matrix B
        for (int32_t j = 0; j < SIZE/2; ++j) {
            MPI_Send(((*sub_b)[j]),SIZE/2, MPI_INT, rank+2, TAG_B, MPI_COMM_WORLD);
            MPI_Recv(((*sub_b)[j+SIZE/2]), SIZE/2, MPI_INT, rank+2, TAG_B, MPI_COMM_WORLD, &stat);
        }
    }
    else{
        for (int32_t j = 0; j < SIZE/2; ++j) {
            MPI_Recv((*sub_b)[j], SIZE/2, MPI_INT, rank-2, TAG_B, MPI_COMM_WORLD, &stat);
            MPI_Send(((*sub_b)[j+SIZE/2]), SIZE/2, MPI_INT, rank-2, TAG_B, MPI_COMM_WORLD);
        }
    }


    printf("Estoy por procesar - rank:%d\n",rank );

    int32_t tmp;
    // PROUCTO dejando fija la columna de B
    for (int32_t i=0; i<SIZE/2; i++){ //i para las filas de la matriz resultante
        for (int32_t j=0; j<SIZE/2 ;j++){ // i para las columnas de la matriz resultante
            tmp = 0 ;
            for (int32_t k=0; k<SIZE; k++){ //k para realizar la multiplicacion de los elementos
                tmp += (*sub_a)[j][k] * (*sub_b)[k][i];
            }
            (*sub_c)[j][i] = tmp;
        }
    }

    if(rank == MASTER) {
        for (int32_t j=0; j<SIZE/2; ++j) {
            MPI_Recv(&((*matrix_c)[j][SIZE/2]), SIZE/2, MPI_INT, 1, TAG_C, MPI_COMM_WORLD, &stat);
            MPI_Recv(&((*matrix_c)[j+SIZE/2]), SIZE/2, MPI_INT, 2, TAG_C, MPI_COMM_WORLD, &stat);
            MPI_Recv(&((*matrix_c)[j+SIZE/2][SIZE/2]), SIZE/2, MPI_INT, 3, TAG_C, MPI_COMM_WORLD, &stat);

//            MPI_Send(sub_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD);
//            MPI_Recv(matrix_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD, &stat);
        }

        if(rank==0){ // el master copia los valores del producto a la matriz original
          for (int32_t i = 0; i < SIZE/2 ; ++i) {
              for (int32_t j = 0; j <SIZE/2 ; j++) {
                (*matrix_c)[i][j]=(*sub_c)[i][j];
              }
          }
        }

        //IMPRIME C
/*        for (int32_t i = 0; i < SIZE; ++i) {
            for (int j = 0; j < SIZE; ++j) {
                printf("%d ", matrix_c[i][j]);
            }
            printf("\n");
        }
*/
        // COMPROBACION
        if(control((int32_t**)matrix_c,SIZE*SIZE)==SUCCESS)
            printf("TOdo OK wacho!\n");
        else
            printf("SE PUDRIO TODO!\n");
    }
    else {
        for (int32_t j = 0; j < SIZE/2; ++j) {
            MPI_Send((*sub_c)[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

int control(int32_t **matrix, int32_t size)
{
  int32_t* aux= (int32_t*)matrix;
  for(int32_t i=0;i<size-1;i++){
    //printf("%d ", *(((int*)matrix)+i));
    if(*(aux+i)>=*(aux+i+1)) return ERROR;
  }
    return SUCCESS;
}

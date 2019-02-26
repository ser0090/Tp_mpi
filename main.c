#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 1024
#define MASTER 0
#define TAG_A 1
#define TAG_B 2
#define TAG_C 3
#define ERROR 1
#define SUCCESS 0


int control(int **matrix, int size);

int main(int argc, char** argv) {

    int world_size;
    int rank;

    int matrix_a[SIZE][SIZE];
    int matrix_b[SIZE][SIZE];
    int matrix_c[SIZE][SIZE];

    int sub_a[SIZE/2][SIZE];
    int sub_c[SIZE/2][SIZE/2];
    int sub_b[SIZE][SIZE/2];
    printf("antes del mpi init\n");
    MPI_Init(NULL, NULL);   // Initialize the MPI environment
    MPI_Status stat;        // required variable for receive routines

    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // Get the rank of the process
    printf("por hacer algo rankl %d\n", rank);
    if(rank == MASTER){
        for (int i = 0; i < SIZE; ++i) {
            for (int j = 0; j < SIZE; ++j) {
                matrix_a[i][j]=j+i*SIZE;
                matrix_b[i][j]=0;
                //printf("%d ",matrix_a[i][j]);
            }
            printf("columna %d",i);
        }
        for (int i = 0; i < SIZE; ++i)
            matrix_b[i][i]=1; //matriz identidad

    }
    printf("Estoy por enviar/recibir primer transmision rank:%d\n",rank );

    for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < 2 ; ++i) {
            for (int j = 0; j <SIZE/2 ; j++) {
                if (rank == MASTER && (i+2*k)){
                    MPI_Send(&matrix_a[j+k*SIZE/2][i*SIZE/2], SIZE/2, MPI_INT, (i+2*k), TAG_A, MPI_COMM_WORLD);
                    MPI_Send(&matrix_b[j+k*SIZE/2][i*SIZE/2], SIZE/2, MPI_INT, (i+2*k), TAG_B, MPI_COMM_WORLD);
                }
                if(rank == (i+2*k) && (i+2*k)){
                    MPI_Recv(&sub_a[j][(rank&1)*SIZE/2], SIZE/2, MPI_INT, MASTER, TAG_A, MPI_COMM_WORLD, &stat);
                    //TODO: revisar j-(rank&1) esta medio hardcode
                    MPI_Recv(&sub_b[j-(rank&1)+k*SIZE/2][i*SIZE/2], SIZE/2, MPI_INT, MASTER, TAG_B, MPI_COMM_WORLD, &stat);
                }
            }
        }
    }

    if(rank==0){ // el master copia los valores a la sub_matrix
      printf("Estoy por copiar mi sub rank:%d\n",rank );
      for (int i = 0; i < SIZE/2 ; ++i) {
          for (int j = 0; j <SIZE/2 ; j++) {
            sub_a[i][j]=matrix_a[i][j];
            sub_b[i][j]=matrix_b[i][j];
          }
      }
    }

    printf("Estoy por enviar sub rank:%d\n",rank );

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


    printf("Estoy por procesar rank:%d\n",rank );

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

//            MPI_Send(sub_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD);
//            MPI_Recv(matrix_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD, &stat);
        }

        if(rank==0){ // el master copia los valores del producto a la matriz original
          for (int i = 0; i < SIZE/2 ; ++i) {
              for (int j = 0; j <SIZE/2 ; j++) {
                matrix_c[i][j]=sub_c[i][j];
              }
          }
        }

        //IMPRIME C
/*        for (int i = 0; i < SIZE; ++i) {
            for (int j = 0; j < SIZE; ++j) {
                printf("%d ", matrix_c[i][j]);
            }
            printf("\n");
        }
*/
        // COMPROBACION
        if(control((int**)matrix_c,SIZE*SIZE)==SUCCESS)
            printf("TOdo OK wacho!\n");
        else
            printf("SE PUDRIO TODO!\n");
    }
    else {
        for (int j = 0; j < SIZE/2; ++j) {
            MPI_Send(sub_c[j], SIZE/2, MPI_INT, MASTER, TAG_C, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

int control(int **matrix, int size)
{
  int* aux= (int*)matrix;
  for(int i=0;i<size-1;i++){
    //printf("%d ", *(((int*)matrix)+i));
    if(*(aux+i)>=*(aux+i+1)) return ERROR;
  }
    return SUCCESS;
}

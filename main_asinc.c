#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/sysinfo.h>

#define SIZE 1700
#define MASTER 0
#define TAG_A 1
#define TAG_B 2
#define TAG_C 3
#define ERROR 1
#define SUCCESS 0
#define SUB_SIZE SIZE/2

typedef int32_t my_matrix_t [SIZE][SIZE];

typedef int32_t my_sub_a_t [SUB_SIZE][SIZE];
typedef int32_t my_sub_b_t [SIZE][SUB_SIZE];
typedef int32_t my_sub_c_t [SUB_SIZE][SUB_SIZE];

int control(int32_t **matrix, int32_t size);

int main(int argc, char** argv) {

    int32_t world_size;
    int32_t rank;

    int32_t procs_count;
    double start_time;
    my_matrix_t* matrix_a;
    my_matrix_t* matrix_b;
    my_matrix_t* matrix_c;

    my_sub_c_t* sub_c;
    my_sub_a_t* sub_a;
    my_sub_b_t* sub_b;

    MPI_Request *fragment_a_req[4]={NULL,NULL,NULL,NULL};
    MPI_Request *fragment_b_req[4]={NULL,NULL,NULL,NULL};
    MPI_Request *fragment_c_req[4]={NULL,NULL,NULL,NULL};


    MPI_Init(NULL, NULL);   // Initialize the MPI environment
    MPI_Status stat;        // required variable for receive routines

    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // Get the rank of the process

/*#################################################################
######################  INICIALIZACION  ###########################
##################################################################*/

    sub_a = (my_sub_a_t*) malloc((SUB_SIZE)*SIZE*sizeof(int32_t));
    sub_b = (my_sub_b_t*) malloc(SIZE*(SUB_SIZE)*sizeof(int32_t));
    sub_c = (my_sub_c_t*) malloc((SUB_SIZE)*(SUB_SIZE)*sizeof(int32_t));

    if(rank == MASTER){
	double start=omp_get_wtime();
	start_time=start;
        matrix_a = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));
        matrix_b = (my_matrix_t*) calloc(SIZE*SIZE, sizeof(int32_t));
        matrix_c = (my_matrix_t*) malloc(SIZE*SIZE*sizeof(int32_t));

        for (int32_t i = 0; i < SIZE; ++i) {
            for (int32_t j = 0; j < SIZE; ++j) {
                (*matrix_a)[i][j]=j+i*SIZE;
                //(*matrix_b)[i][j]=0;
                //printf("%d ",matrix_a[i][j]);
            }
        }
        for (int32_t i = 0; i < SIZE; ++i) {
            (*matrix_b)[i][i] = 1; //matriz identidad
        }
        // el master copia los valores a la sub_matrix
        for (int32_t i = 0; i < SUB_SIZE ; ++i) {
            for (int32_t j = 0; j <SUB_SIZE ; j++) {
                (*sub_a)[i][j] = (*matrix_a)[i][j];
                (*sub_b)[i][j] = (*matrix_b)[i][j];
            }
        }
	printf("tiempo para iniciar la matriz: %f\n",omp_get_wtime()-start);
    }

/*#################################################################
#################  FRACCIONAMIENTO DE MATRICES  ###################
##################################################################*/
    if(rank==MASTER) {
      //fragment_a_req[0]=(MPI_Request*)calloc(SIZE,sizeof(MPI_Request));
/*      fragment_a_req[1]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      fragment_a_req[2]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      fragment_a_req[3]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      //fragment_b_req[0]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      fragment_b_req[1]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      fragment_b_req[2]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      fragment_b_req[3]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));

      //fragment_a_stat[0]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      fragment_a_stat[1]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      fragment_a_stat[2]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      fragment_a_stat[3]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      //fragment_b_stat[0]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      fragment_b_stat[1]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      fragment_b_stat[2]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      fragment_b_stat[3]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
*/    MPI_Request aux_req;
       double start = omp_get_wtime();
        for (int32_t k = 0; k < 2; ++k) {
            for (int32_t i = 0; i < 2; ++i) {
                if ((i + 2 * k)) {
                    for (int32_t j = 0; j < SIZE / 2; j++) {
  //                    printf("por enviar a rank %d\n", (i+2*k));
                        MPI_Isend(&((*matrix_a)[j + k * SIZE / 2][i * SIZE / 2]), SIZE / 2, MPI_INT, (i + 2 * k), TAG_A,
                                 MPI_COMM_WORLD,&aux_req );//&(fragment_a_req[(i + 2 * k)][j]));
                        MPI_Request_free(&aux_req);
                        MPI_Isend(&((*matrix_b)[j + k * SIZE / 2][i * SIZE / 2]), SIZE / 2, MPI_INT, (i + 2 * k), TAG_B,
                                 MPI_COMM_WORLD, &aux_req );// &(fragment_b_req[(i + 2 * k)][j]));
                        MPI_Request_free(&aux_req);

  //                    printf("enviado a rank %d\n", (i+2*k));
                    }
                }
            }
        }


      //no es necesario que espere por ninguno. cada hilo espera independientemente a tener sus datos antes de continuar.
    //MPI_Waitall(SUB_SIZE,fragment_a_req[1],fragment_a_stat[1]);
    //MPI_Waitall(SUB_SIZE,fragment_a_req[2],fragment_a_stat[2]);
    //MPI_Waitall(SUB_SIZE,fragment_a_req[3],fragment_a_stat[3]);
    //MPI_Waitall(SUB_SIZE,fragment_b_req[1],fragment_b_stat[1]);
    //MPI_Waitall(SUB_SIZE,fragment_b_req[2],fragment_b_stat[2]);
    //MPI_Waitall(SUB_SIZE,fragment_b_req[3],fragment_b_stat[3]);

    printf("tiempo en fraccionar: %f\n",omp_get_wtime()-start);

    }
    else{
      fragment_a_req[rank]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      fragment_b_req[rank]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
//      fragment_a_stat[rank]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
//      fragment_b_stat[rank]=(MPI_Status*)malloc(SUB_SIZE*sizeof(MPI_Status));
      if(fragment_a_req[rank]==NULL || fragment_b_req[rank]==NULL){
        printf("CALLOC_ERROR\n"); return 1;
      }

      for (int32_t k = 0; k < 2; ++k) {
          for (int32_t i = 0; i < 2 ; ++i) {
              if(rank == (i+2*k) && (i+2*k)){ //rank == (i+2*k) es para que cada proceso solo cargue su cuadrante, sino los cargaria a todos.
                  for (int32_t j = 0; j <SUB_SIZE ; j++) {
                      MPI_Irecv(&((*sub_a)[j][(rank&1)*SUB_SIZE]), SUB_SIZE, MPI_INT, MASTER, TAG_A, MPI_COMM_WORLD, &(fragment_a_req[rank][j]));
                      //TODO: revisar j-(rank&1) esta medio hardcode
                      MPI_Irecv(&((*sub_b)[j-(rank&1)+k*SUB_SIZE][i*SUB_SIZE]), SUB_SIZE, MPI_INT, MASTER, TAG_B, MPI_COMM_WORLD, &(fragment_b_req[rank][j]));
                  }
              }
          }
      }
      printf("esperando transferencia - rank %d\n",rank);

      MPI_Waitall(SUB_SIZE,fragment_a_req[rank],MPI_STATUS_IGNORE);
      MPI_Waitall(SUB_SIZE,fragment_b_req[rank],MPI_STATUS_IGNORE);

    }

//    printf("Estoy por enviar sub rank:%d\n",rank );

    if(!(rank/2)){ // Pasaje de Matrix B
      fragment_b_req[rank+2]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
        MPI_Request aux_req;
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            MPI_Isend(((*sub_b)[j]),SUB_SIZE, MPI_INT, rank+2, TAG_B, MPI_COMM_WORLD,&aux_req);
            MPI_Request_free(&aux_req);
            MPI_Irecv(((*sub_b)[j+SUB_SIZE]), SUB_SIZE, MPI_INT, rank+2, TAG_B, MPI_COMM_WORLD, &(fragment_b_req[rank+2][j]));
        }
    }
    else{
        fragment_b_req[rank-2]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
        MPI_Request aux_req;
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            MPI_Irecv((*sub_b)[j], SUB_SIZE, MPI_INT, rank-2, TAG_B, MPI_COMM_WORLD, &(fragment_b_req[rank-2][j]));
            MPI_Isend(((*sub_b)[j+SUB_SIZE]), SUB_SIZE, MPI_INT, rank-2, TAG_B, MPI_COMM_WORLD,&aux_req);
            MPI_Request_free(&aux_req);
        }
    }

    if(!(rank&1)){ //is rank es par para A
        fragment_a_req[rank+1]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
        MPI_Request aux_req;
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            MPI_Isend(&((*sub_a)[j]),SUB_SIZE, MPI_INT, rank+1, TAG_A, MPI_COMM_WORLD,&aux_req);
            MPI_Request_free(&aux_req);
            MPI_Irecv(&((*sub_a)[j][SUB_SIZE]), SUB_SIZE, MPI_INT, rank+1, TAG_A, MPI_COMM_WORLD, &(fragment_a_req[rank+1][j]));
        }
    }
    else{
      fragment_a_req[rank-1]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
      MPI_Request aux_req;
        for (int32_t j = 0; j < SUB_SIZE; ++j) {
            MPI_Irecv(&((*sub_a)[j]), SUB_SIZE, MPI_INT, rank-1, TAG_A, MPI_COMM_WORLD, &(fragment_a_req[rank-1][j]));
            MPI_Isend(&((*sub_a)[j][SUB_SIZE]), SUB_SIZE, MPI_INT, rank-1, TAG_A, MPI_COMM_WORLD,&aux_req);
            MPI_Request_free(&aux_req);
        }
    }





/*#################################################################
#######################   PROCESAMIENTO   #########################
##################################################################*/
    //printf("Estoy por procesar - rank:%d\n",rank );

    register int32_t tmp;
    procs_count = get_nprocs();
//    omp_set_num_threads(procs_count/2);

    // PROUCTO dejando fija la columna de B
   // #pragma openmp parallel
//{
//   #pragma omp parallel for
    int32_t wait_b_from = (!(rank/2))? rank+2:rank-2;
    int32_t wait_a_from = (!(rank&1))? rank+1:rank-1;
    MPI_Waitall(SUB_SIZE,fragment_b_req[wait_b_from],MPI_STATUS_IGNORE);

    //despliegue del procesamiento para esperar lo datos de a
    for (int32_t j=0; j<SUB_SIZE ;j++){ // i para las columnas de la matriz resultante
        MPI_Wait(&(fragment_a_req[wait_a_from][j]),MPI_STATUS_IGNORE);
        tmp = 0 ;
        for (int32_t k=0; k<SIZE; k++){ //k para realizar la multiplicacion de los elementos
            tmp += (*sub_a)[j][k] * (*sub_b)[k][0];
        }
        (*sub_c)[j][0] = tmp;
    }


    for (int32_t i=1; i<SUB_SIZE; i++){ //i para las filas de la matriz resultante
        for (int32_t j=0; j<SUB_SIZE ;j++){ // i para las columnas de la matriz resultante
            tmp = 0 ;
            for (int32_t k=0; k<SIZE; k++){ //k para realizar la multiplicacion de los elementos
                tmp += (*sub_a)[j][k] * (*sub_b)[k][i];
            }
            (*sub_c)[j][i] = tmp;
        }
    }
//}

/*#################################################################
#######################   RECONSTRUCCION   #########################
##################################################################*/
    if(rank == MASTER) {
        for(int i=1;i<4;i++){
              fragment_c_req[i]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
          }
	      double start = omp_get_wtime();
        for (int32_t j=0; j<SUB_SIZE; ++j) {
            MPI_Irecv(&((*matrix_c)[j][SUB_SIZE]), SUB_SIZE, MPI_INT, 1, 1, MPI_COMM_WORLD, &(fragment_c_req[1][j]));
            MPI_Irecv(&((*matrix_c)[j+SUB_SIZE]), SUB_SIZE, MPI_INT, 2, 2, MPI_COMM_WORLD, &(fragment_c_req[2][j]));
            MPI_Irecv(&((*matrix_c)[j+SUB_SIZE][SUB_SIZE]), SUB_SIZE, MPI_INT, 3, 3, MPI_COMM_WORLD, &(fragment_c_req[3][j]));

        }

        for (int32_t i = 0; i < SUB_SIZE ; ++i) {
            for (int32_t j = 0; j <SUB_SIZE ; j++) {
              (*matrix_c)[i][j]=(*sub_c)[i][j];
            }
        }
        printf("TERMINO TODO - rank %d\n",rank);
        MPI_Waitall(SUB_SIZE,fragment_c_req[1],MPI_STATUS_IGNORE);
        MPI_Waitall(SUB_SIZE,fragment_c_req[2],MPI_STATUS_IGNORE);
        MPI_Waitall(SUB_SIZE,fragment_c_req[3],MPI_STATUS_IGNORE);

        printf("tiempo en construir c: %f\n",omp_get_wtime()-start);
        printf("tiempo Total: %f\n",omp_get_wtime()-start_time);
        // COMPROBACION
        if(control((int32_t**)matrix_c,SIZE*SIZE)==SUCCESS)
            printf("Todo OK wacho ^^\n");
        else
            printf("SE PUDRIO TODO T.T\n");

        free(matrix_a);
        free(matrix_b);
        free(matrix_c);
    }
    else {
        printf("termino de procesar - rank %d\n",rank);
        //MPI_Request aux_req;
        fragment_c_req[rank]=(MPI_Request*)malloc(SUB_SIZE*sizeof(MPI_Request));
        for (int32_t j = 0; j < SUB_SIZE ; ++j) {
            //printf("enviando fila %d - rank %d\n",j,rank);
            MPI_Isend((*sub_c)[j], SUB_SIZE, MPI_INT, MASTER, rank, MPI_COMM_WORLD,&(fragment_c_req[rank][j]));
            //printf("sent fila %d - rank %d\n",j,rank);
            //MPI_Request_free(&aux_req);
        }
        printf("TERMINO TODO - rank %d\n",rank);

        MPI_Waitall(SUB_SIZE,fragment_c_req[rank],MPI_STATUS_IGNORE);
    }
/*
    free(sub_a);
    free(sub_b);
    free(sub_c);

    for(int i=0;i<4;i++){
      free(fragment_a_req[i]);
      free(fragment_b_req[i]);
      free(fragment_a_stat[i]);
      free(fragment_b_stat[i]);
    }

*/
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

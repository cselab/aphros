#include <mpi.h>
#include "ParIO.h" 

int main(int argc, char *argv[])
{
	int rank, size;

	MPI_Init(&argc,&argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size); 

	float data[32];

	for (int i = 0; i < 32; i++)
		data[i] = rank*1000 + i;

	MPI_ParIO pio;


	pio.Init((char *)"output.bin", 32*4);

	int id = (rank + 1)%size;

	pio.Write((char *)data, id);

	pio.Finalize();

	MPI_Finalize();
	
	return 0;
} 

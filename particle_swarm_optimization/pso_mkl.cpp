#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<mpi.h>

#define PI 3.141592653589
#define pnum 200
#define dim 1024
#define iternum 500000
#define low -5.12
#define high 5.12
#define vmax 2000
#define c1 2
#define c2 2
#define w 0.5
#define alpha 1


double Rastrigin(double a[]) {
	double sum = 0.0;
	for (int i = 0;i < dim;i++)
	{
		sum += a[i] * a[i] - 8 * cos(2 * PI*a[i]) + 8;
	}
	return sum;
}

double splras(unsigned iter, int my_rank) {
	double gfit, particle[pnum][dim], local_best[pnum][dim], local_fit[pnum], global_best[dim], particle_v[pnum][dim], particle_fit[pnum];
	int i, j;
	srand(pow(my_rank + 1, my_rank + 1));
	for (i = 0; i<pnum; i++) {
		for (j = 0; j<dim; j++) {
			particle[i][j] = low + (high - low)*1.0*rand() / RAND_MAX;
			local_best[i][j] = particle[i][j];
			particle_v[i][j] = -vmax + 2 * vmax*1.0*rand() / RAND_MAX;
		}
	}
	for (i = 0; i<pnum; i++) {
		particle_fit[i] = Rastrigin(particle[i]);
		local_fit[i] = particle_fit[i];
	}
	gfit = local_fit[0];
	j = 0;
	for (i = 1; i<pnum; i++) {
		if (local_fit[i]<gfit) {
			gfit = local_fit[i];
			j = i;
		}
	}
	for (i = 0; i<dim; i++) {
		global_best[i] = local_best[j][i];
	}
	for (int k = 0;i < iter;i++) {
		for (i = 0; i<pnum; i++) {
			for (j = 0; j<dim; j++) {
				particle[i][j] += alpha*particle_v[i][j];
				if (particle[i][j] > high) {
					particle[i][j] = high;
				}
				if (particle[i][j] < low) {
					particle[i][j] = low;
				}
			}
		}
		for (i = 0; i<pnum; i++) {
			particle_fit[i] = Rastrigin(particle[i]);
			if (particle_fit[i] < local_fit[i]) {
				local_fit[i] = particle_fit[i];
				for (j = 0; j<dim; j++) {
					local_best[i][j] = particle[i][j];
				}
			}
		}
		for (i = 0, j = -1; i<pnum; i++) {
			if (local_fit[i]<gfit) {
				gfit = local_fit[i];
				j = i;
			}
		}
		if (j != -1) {
			for (i = 0; i<dim; i++) {
				global_best[i] = local_best[j][i];
			}
		}
		for (i = 0; i<pnum; i++) {
			for (j = 0; j<dim; j++) {
				particle_v[i][j] = w*particle_v[i][j] +
					c1*1.0*rand() / RAND_MAX*(local_best[i][j] - particle[i][j]) +
					c2*1.0*rand() / RAND_MAX*(global_best[j] - particle[i][j]);
				if (particle_v[i][j] > vmax) {
					particle_v[i][j] = vmax;
				}
				if (particle_v[i][j] < -vmax) {
					particle_v[i][j] = -vmax;
				}
			}
		}
	}
	return gfit;
}

int main(int argc, char* argv[])
{
	clock_t start = clock();
	int my_rank, p;
	double start_time, end_time, elapse_time, as;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	unsigned iter = iternum / p;
	start_time = MPI_Wtime();
	as = splras(iter,my_rank);
	end_time = MPI_Wtime();
	elapse_time = end_time - start_time;
	printf("Time on processor%d is %f seconds.\n", my_rank, elapse_time);
	printf("Minimum value on processor%d is %f.\n", my_rank, as);
	MPI_Finalize();
	clock_t stop = clock();
	printf("The total time is %f seconds.\n", (float)(stop - start) / CLOCKS_PER_SEC);
	return 0;
}
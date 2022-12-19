#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include "functions.h"


void
reduce_sum(int p, double *a = nullptr, int n = 0);
void reduce_sum(int p, double *a, int n)
{
  static pthread_mutex_t m     = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t  c_in  = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t  c_out = PTHREAD_COND_INITIALIZER;
  static int t_in  = 0;
  static int t_out = 0;
  static double *r = nullptr;
  int i;

  if (p <= 1)
    return;
  pthread_mutex_lock(&m);
  if (r == nullptr)
    r = a;
  else
    for (i = 0; i < n; i++)
      r[i] += a[i];
  t_in++;
  if (t_in >= p)
  {
    t_out = 0;
    pthread_cond_broadcast(&c_in); // оповещение что все потоки вошли
  }
  else
    while (t_in < p)
      pthread_cond_wait(&c_in, &m);
  if (r != a)
    for (i = 0; i < n; i++)
      a[i] = r[i];
  t_out++;
  if (t_out >= p)
  {
    t_in = 0;
    r = nullptr;
    pthread_cond_broadcast(&c_out);
  }
  else
    while (t_out < p)
      pthread_cond_wait(&c_out, &m);
  pthread_mutex_unlock(&m);
}

double matrix_norm(double *matrix, int n)
{
  double norm = 0, tmp_sum = 0;
  for (int i = 0; i < n; i++)
  {
    tmp_sum = 0;
    for (int j = 0; j < n; j++)
      tmp_sum += fabs(matrix[i * n + j]);
    if (tmp_sum > norm)
      norm = tmp_sum;
  }
  return norm;
}

double Nevyazka(double *A, double *x, double *b, int n)
{
    double *Axb;
    Axb = new double [n];
    for (int i=0; i<n; i++)
        Axb[i]=0;
    double n1 = 0, n2 = 0;

    for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    {
                    Axb[i] += A[i * n + j] * x[j];

                    }
                Axb[i] -= b[i];
                n1 += pow(Axb[i],2);

                n2 += pow(b[i],2);
            }

    n1 = sqrt(n1);
    n2 = sqrt(n2);
    delete[] Axb;
    return n1/n2;
}

double sin_norm(double *x, int n)
{
    double sum = 0;
    int t;
    for(int i = 0; i < n; i++)
    {
        if (i % 2 == 0)
            t = 1;
        else
            t = 0;
        sum += pow (x[i] - t, 2);
    }

    return sqrt(sum);
}

//-------------------------------------------------------------------------

struct thread_data{
    double* A;
    double* b;
    double* Max_ara;
    double* Max_row;
    double* Max_col;
    double* maxi;
	int* indicies;
	int n;
	int index;
	sem_t* sem_1;
	sem_t* sem_2;
	pthread_mutex_t* mutex;
	int* fl;
	int* thread_finished_counter;
	int number_of_threads;
};

int Jora_parallel(int n, double* A, double* b, double* x, int* indicies, int number_of_threads){

    pthread_t *threads;
    struct thread_data *threads_data;

    sem_t sem_1;
	sem_init(&sem_1, 0, 0);

	sem_t sem_2;
	sem_init(&sem_2, 0, 0);

    double* Max_ara;
    double* Max_col;
    double* Max_row;

    Max_ara = new double [number_of_threads];
    Max_row = new double [number_of_threads];
    Max_col = new double [number_of_threads];

    double maxi = -10000000;
    int fl = 1;
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, NULL);

       threads = new pthread_t[number_of_threads];
       // pthread_t threads[number_of_threads];
	int thread_finished_counter = 0;

        threads_data = new struct thread_data[number_of_threads];
      //  struct thread_data threads_data[number_of_threads];
	for (int i = 0; i < number_of_threads; i++){
		threads_data[i].A = A;
		threads_data[i].Max_ara = Max_ara;
		threads_data[i].Max_row = Max_row;
		threads_data[i].Max_col = Max_col;
		threads_data[i].b = b;
        threads_data[i].maxi = &maxi;
		threads_data[i].n = n;
        threads_data[i].indicies = indicies;
		threads_data[i].sem_1 = &sem_1;
		threads_data[i].fl = &fl;
		threads_data[i].sem_2 = &sem_2;
		threads_data[i].mutex = &mutex;
		threads_data[i].thread_finished_counter = &thread_finished_counter;
		threads_data[i].index = i;
		threads_data[i].number_of_threads = number_of_threads;
	}

	for (int i = 1; i < number_of_threads; i++){
                pthread_create(&threads[i], NULL, Jora_thread_func,(void*) &threads_data[i]);
                //std::cout<<"\nN - "<< thread_finished_counter<<std::endl;

	}
    Jora_thread_func(threads_data + 0);
	for (int i = 1; i < number_of_threads; i++){

		pthread_join(threads[i], NULL);
        //std::cout<<"\nM - "<< thread_finished_counter<<std::endl;
	}


    if (fl == -69)

        return -1;

    for (int i = 0; i < n; i++)
		x[indicies[i]] = b[i];
    delete [] threads;
    delete [] threads_data;

return 0;
}

void* Jora_thread_func(void* arg){
	struct thread_data* data = (struct thread_data*)arg;
    int max_row, max_col, tr;
    double tmp;
    double max_el;
    double* A = data->A;
    double* b = data->b;


   // double max_el = *data->maxi;
	int* indicies = data->indicies;
	int iloveprogat;
	int n = data->n;
   	int p = data->index;
    int N = data->number_of_threads;
    double norm = matrix_norm(A, n);

    double *Max_ara = data->Max_ara;
    double *Max_row = data->Max_row;
    double *Max_col = data->Max_col;



      int i_min = p * (n / N);
		int i_max;

        if (p == N - 1)
            i_max = n;
        else
            i_max = (p + 1) * (n / N);





	for (int i = 0; i < n; i++) {
        for (int g = 0; g < N; g++)
             Max_ara[g] = fabs(A[i * n + i]);
         //printf("%d",i);
		/*pthread_mutex_lock(data->mutex);

		(*data->thread_finished_counter)++;

		if(*data->thread_finished_counter == data->number_of_threads)
        {

            max_el = fabs(A[i * n + i]);
            max_row = i;
            max_col = i;

            for (int j = i; j < n; j++)
                for (int k = i; k < n; k++)
                    if (max_el < fabs(A[j * n + k]))
                    {
                        max_el = fabs(A[j * n + k]);
                        max_row = j;
                        max_col = k;
                    }
            /*if (fabs(A[max_row * n + max_col]) < 1.2e-16 * norm)
                std::cout<<max_row<<std::endl;
            if (fabs(A[max_row * n + max_col]) <= 1.2e-16 * norm){

                *data->fl = -69;
                for (int h = 0; h < data->number_of_threads - 1; h++){
                    sem_post(data->sem_1);
                    sem_wait(data->sem_2);
                  }
                pthread_mutex_unlock(data->mutex);

                return NULL;}

            tr = indicies[i];
            indicies[i] = indicies[max_col];
            indicies[max_col] = tr;

            tmp = b[i];
            b[i] = b[max_row];
            b[max_row] = tmp;

            for (int j = 0; j < n; j++)
            {
                tmp = A[i * n + j];
                A[i * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }

            for (int j = 0; j < n; j++)
            {
                tmp = A[j * n + i];
                A[j * n + i] = A[j * n + max_col];
                A[j * n + max_col] = tmp;
            }

            tmp = 1.0 / A[i * n + i];
            for (int j = i; j < n; j++)
                A[i * n + j] *= tmp;
            b[i] *= tmp;

            //std::cout<<std::endl;
           // print_matrix(n, n, n, A);
            //std::cout<<std::endl;


            for (int i = 0; i < data->number_of_threads - 1; i++){
                sem_post(data->sem_1);
                //sem_wait(data->sem_2);
            }
            *data->thread_finished_counter = 0;
        }
		else {
			pthread_mutex_unlock(data->mutex);
			sem_wait(data->sem_1);
			//sem_post(data->sem_2);
			pthread_mutex_lock(data->mutex);
		}

		pthread_mutex_unlock(data->mutex);*/


        //max_el = fabs(A[i * n + i]);
        //max_row = i;
       // max_col = i;

        int i_start = i + p * ((n-i)/N);
        int i_end = i + (p+1) * ((n-i)/N);
        if (p == N - 1)
            i_end = n;


		/*int a = i + 1;
		int b = sle->size;
		int p = data->index;
		int N = data->number_of_threads;

		int i_min = a + p * ((b - a) / N);
		int i_max = p == N - 1 ? b : a + (p + 1) * ((b - a) / N);*/


        for (int j = i_start; j < i_end; j++)
            for (int k = i; k < n; k++)
                if (Max_ara[p] < fabs(A[j * n + k]))
                {
                    Max_ara[p] = fabs(A[j * n + k]);
                    Max_row[p] = j;
                    Max_col[p] = k;
                }




        reduce_sum(N);


        max_el = Max_ara[0];
        for (int y = 0; y < p; y++)
        {
            if (Max_ara[y] > max_el)
                {max_el = Max_ara[y];
                iloveprogat = y;}
        }

        max_row = Max_row[iloveprogat];
        max_col = Max_col[iloveprogat];

        if (p == 0)
        {


            if (fabs(A[max_row * n + max_col]) <= 1.2e-16 * norm){

                *data->fl = -69;
                 reduce_sum(N);
                 //for (int h = 0; h < data->number_of_threads - 1; h++){

                    //sem_post(data->sem_1);
                    //sem_wait(data->sem_2);
                  //}
               // pthread_mutex_unlock(data->mutex);

                return NULL;}

            tr = indicies[i];
            indicies[i] = indicies[max_col];
            indicies[max_col] = tr;

            tmp = b[i];
            b[i] = b[max_row];
            b[max_row] = tmp;

            for (int j = 0; j < n; j++)
            {
                tmp = A[i * n + j];
                A[i * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }

            for (int j = 0; j < n; j++)
            {
                tmp = A[j * n + i];
                A[j * n + i] = A[j * n + max_col];
                A[j * n + max_col] = tmp;
            }

            tmp = 1.0 / A[i * n + i];
            for (int j = i; j < n; j++)
                A[i * n + j] *= tmp;
            b[i] *= tmp;
        }



        reduce_sum(N);


        if (*data->fl == - 69)
            return NULL;





		for (int j = i_min; j < i_max; j++)
        {


			if (j<i)
			{tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

            if (j>i){
			tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

			else{}

        }

	}




	return NULL;
}




int Jordan_mid(int n, double* A, double* b, double* x, int* index)
{
	int max_row, max_col, tr;
	double tmp, max_el;

    double norm = matrix_norm(A, n);

	for (int i = 0; i < n; i++)
		index[i] = i;

	for (int i = 0; i < n; i++)
	{
		max_el = fabs(A[i * n + i]);
		max_row = i;
		max_col = i;

		for ( int j = i; j < n; j++)
			for (int k = i; k < n; k++)
				if (max_el < fabs(A[j * n + k]))
				{
					max_el = fabs(A[j * n + k]);
					max_row = j;
					max_col = k;
				}

		if (fabs(A[max_row * n + max_col]) < 1.2e-16 * norm)
			return -1;

		tr = index[i];
		index[i] = index[max_col];
		index[max_col] = tr;

		tmp = b[i];
		b[i] = b[max_row];
		b[max_row] = tmp;

		for (int j = 0; j < n; j++)
		{
			tmp = A[i * n + j];
			A[i * n + j] = A[max_row * n + j];
			A[max_row * n + j] = tmp;
		}

        for (int j = 0; j < n; j++)
		{
			tmp = A[j * n + i];
			A[j * n + i] = A[j * n + max_col];
			A[j * n + max_col] = tmp;
		}

		tmp = 1.0 / A[i * n + i];
		for (int j = i; j < n; j++)
			A[i * n + j] *= tmp;
		b[i] *= tmp;

		for (int j = 0; j < n; j++)
		{
			if (j<i)
			{tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

            if (j>i){
			tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

			else{}
		}

                //std::cout<<"\n--------------------\n";
                //print_matrix(n,n,n, A);

        }



	for (int i = 0; i < n; i++)
		x[index[i]] = b[i];

	return 0;}

/*int main(void)
{
double a[]= {3, 2, -5, 2, -1, 3, 1, 2, -1};
double b[] = {-1, 13 , 9};
double x[] = {0, 0, 0};
int index[] = {0, 0, 0, 0, 0, 0};
int ara;
ara = Jordan_mid(3, a, b, x ,index);
std::cout<<ara<<"\n";
std::cout<<x[0]<<"\n";
std::cout<<x[1]<<"\n";
std::cout<<x[2]<<"\n";
return 1;
}*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void print_matrix(int n, int m, double **matrix){
	for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			printf("%lf ",matrix[i][j]);
		}
		
		printf("\n");
	}
}


double* thomas_algo(int n, double **matrix, double *d){
	int i=0;
	int temp_n=n;
	
	double *w = (double*)malloc((n-1)*sizeof(double));
	double *g = (double*)malloc(n*sizeof(double));
	
	w[0] = matrix[0][1]/matrix[0][0];
	g[0] = d[0]/matrix[0][0];
	
	for (i=1;i<n-1;i++){
		w[i] = matrix[i][i+1]/(matrix[i][i] - matrix[i][i-1] * w[i-1]);
	}
	
	for (i=1;i<n;i++){
		g[i] = (d[i] - matrix[i][i-1]*g[i-1])/(matrix[i][i] - matrix[i][i-1]*w[i-1]);
	}
		
	d[n-1]=g[n-1];
	
	for (int i=n;i>0;i--){
        d[i-1] = g[i-1] - w[i-1]*d[i];
	}
	
	free(w);
	free(g);
	return d;
}

int main(int argc, char* argv[]){
	int n=0;
	int mode_read=0;
	int mode_check=0;
	int mode_print=0;
	int mode_help=0;
	
	double temp=0;
	double **matrix=NULL;
	double *d=NULL;
	double *tmp_d=NULL;
	
	const char* fin_name=NULL;
	const char* fout_name=NULL;
	
	FILE *fin=NULL;
	FILE *fout=NULL;
	
	/*-- Чтение и разбор параметров --*/
	
	if (argc>2){
		if (strstr(argv[1],".txt")){
			fin_name=argv[1];
		}
		
		if (strstr(argv[2],".txt")){
			fout_name=argv[2];
		}
	} else {
		printf("Error format. Format executable_file input.txt output.txt [options]\n");
		return -1;
	}
	
	for (int i=3;i<argc;i++){
		if (strcmp(argv[i],"-diag")==0){ //форма входных данных
			mode_read=1;
		}
		
		if (strcmp(argv[i],"-check")==0){ //проверить результат после выполнения
			mode_check=1;
		}
		
		if (strcmp(argv[i],"-v")==0){ //выводить на экран
			mode_print=1;
		}
		
		if (strcmp(argv[i],"-h")==0){ //выводить на экран
			mode_help=1;
		}
	}
	
	if (mode_help==1){
		printf("executable_file input.txt output.txt [options]\n");
		printf(" options: \n -diag read from file where only three diagonals of matrix (each diagonal size is n)\n");
		printf(" -check checking result\n");
		printf(" -v print datas on display\n");
		printf(" -h help\n");
		return 1;
	}
	
	
	fin = fopen(fin_name,"r");
	fscanf(fin,"%d",&n);
	
	if (n<3){
		if (mode_print==1){
			printf("n cannot be less than 3. Error. Exit");
		}
		return -1;
	}
	
	/*-- Выделение памяти и чтение данных --*/
	
	matrix=(double**)malloc(n*sizeof(double *));
	d = (double*)malloc(n*sizeof(double));
	tmp_d = (double*)malloc(n*sizeof(double));
	
	for (int i=0;i<n;i++){
		matrix[i]=(double*)malloc(n*sizeof(double));
	}
	
	for (int i=0; i<n; i++){
		for (int j=0; j<n;j++){
			matrix[i][j] = 0;
		}
	}
	
	
	if (mode_read==0){
		for (int i=0; i<n; i++){
			for (int j=0; j<n;j++){
				fscanf(fin,"%lf",&temp);
				matrix[i][j] = temp;
			}
		}
		
		for (int i=0; i<n; i++){
			fscanf(fin,"%lf",&temp);
			d[i]=temp;
			if (mode_check==1){
				tmp_d[i]=temp;
			}
		}

		fclose(fin);
	}
	
	if (mode_read==1){
		/* Считывание данных
		первая строка - под главной диагональю;
		вторая строка - главная диагональю
		третья строка - над главной диагональю
		размер каждой n
		*/
		
		for (int i=0; i<n-1; i++){
			fscanf(fin,"%lf",&temp);
			matrix[i+1][i]=temp;
		}
	
		fscanf(fin,"%lf",&temp);
	
		for (int i=0; i<n; i++){
			fscanf(fin,"%lf",&temp);
			matrix[i][i]=temp;
		}
	
		for (int i=0; i<n-1; i++){
			fscanf(fin,"%lf",&temp);
			matrix[i][i+1]=temp;
		}
	
		fscanf(fin,"%lf",&temp);
	
		for (int i=0; i<n; i++){
			fscanf(fin,"%lf",&temp);
			d[i]=temp;
			if (mode_check==1){
				tmp_d[i]=temp;
			}
		}

		fclose(fin);
	}
	
	if (mode_print==1){
		printf("Ax=B, where A triangular matrix:\n");
		printf("A:\n");
		print_matrix(n,n, matrix);
		printf("B:\n");
		for (int i=0;i<n;i++){
			printf("b%d:%lf\n",i+1,d[i]);
		}
	}
	
	d=thomas_algo(n,matrix,d);
	
	if (mode_print==1){
		printf("X:\n");
		for (int i=0;i<n;i++){
			printf("x%d:%lf\n",i+1,d[i]);
		}
	}
	
	if (mode_check==1){
		double eps=1e-5;
		printf("Check: AX?=B\n");
		double **test_matrix=(double**)malloc(n*sizeof(double *));
	
		for (int i=0;i<n;i++){
			test_matrix[i]=(double*)malloc(n*sizeof(double));
		}
		
		
		for (int i = 0; i < n; i++){
			test_matrix[i][0] = 0;
			for (int k = 0; k < n; k++){
				test_matrix[i][0] += matrix[i][k] * d[k];
			}
		}

		//print_matrix(n,1,test_matrix);
		for (int i=0;i<n;i++){
			if (fabs(test_matrix[i][0]-tmp_d[i])<eps){
				printf("x%d:ok\n",i+1);
			} else {
				printf("x%d:failed\n",i+1);
			}
		}
		
	}
	
	//Вывод в файл
	fout=fopen(fout_name,"w");
	
	fprintf(fout,"%d\n",n);
	for (int i=0;i<n;i++){
		fprintf(fout,"%lf\n",d[i]);
	}
	
	fclose(fout);
	free(d);
	
	for (int i=0;i<n;i++){
		free(matrix[i]);
	}
	
	free(matrix);
	return 0;
}
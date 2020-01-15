
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <malloc.h>

#define MATRIX float *
#define KDTREE struct kdtree_node * // modificare con il tipo di dato utilizzato

void *get_block(int size, int elements)
{
    return _mm_malloc(elements * size, 16);
}

void free_block(void *p)
{
    _mm_free(p);
}

//float euclidean_distance(MATRIX qs, int id_qs, MATRIX ds, int id_ds, int k);

/*

int indexList = 0;



float euclidean_distance(MATRIX qs, int id_qs, MATRIX ds, int id_ds, int k)
{
    float somma = 0;
    int i=0;
    for (; i < k; i++)
    {
        somma = somma + powf(qs[id_qs * k + i] - ds[id_ds * k + i],2);
    }
    return sqrtf(somma);
}

void range_query(params *input)
{
    printf("\n inizio QUERY");

    input->QA = (int *)malloc(input->n * input->nq * sizeof(int));
    int num=0;
    if (input->QA == NULL)
    {
        printf("\nNO MEMORIA");
        exit(1);
    }
    int i;
    for (int id_qs= 0; id_qs <  input->nq; id_qs++){
         for (int id_ds= 0; id_ds <  input->n; id_ds++){
        if (euclidean_distance(input->qs, id_qs,input->ds,id_ds,input->k) <= input->r){
        input->QA[id_qs * input->n + indexList] = id_ds;
        indexList++;
        num=num+1;}
        }
        input->QA[id_qs * input->n + indexList] = -1;
        indexList = 0;
    }
    input->nQA=num;
}
*/

void crea_file(char *filename, void *X, int k, int n)
{
    FILE *fp;
    fp = fopen(filename, "wb");
    if (X != NULL)
    {
        fwrite(&k, 4, 1, fp);
        fwrite(&n, 4, 1, fp);
        for (int i = 0; i < n; i++)
        {
            fwrite(X, 4, k, fp);
            X += 4 * k;
        }
    }
    fclose(fp);
}

/*
MATRIX leggi(char *filename)
{
    FILE *fp;
    int rows, cols, status, i;

    fp = fopen(filename, "rb");

    if (fp == NULL)
    {
        printf("'%s': bad data file name!\n", filename);
        exit(0);
    }

    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);

    float *data = (float *)get_block(sizeof(float), rows * cols);
    status = fread(data, sizeof(float), rows * cols, fp);
    fclose(fp);
    for (int i = 0; i < (rows); i=i+1)
    { 
        printf("\n\n");
        for (int j = 0; j < (cols); i=i+1){
        printf("%f ; ", data[i]);
    }
    }
    free_block(data);
}
*/

//2147483647/8 =268435455,875 numero di byte della memoria
//

int main()
{
    // int max_value=4294967295;  /////4294967295;
    // int max_value2=1094967294;  /////4294967295;

    // int * ppp = (int *)get_block(sizeof(int), max_value2);
    // if (ppp == NULL)
    // {
    //     printf("\n problema grande ");
    //     exit(1);
    // }

    char *filename = "prova1";
    char fname[256];
    int k = 128;
    int n_point = 80003; //(int max=2147483647)
    int n_query = 4003;
    float *ds = (float *)get_block(sizeof(float), n_point * k);
    float *qs = (float *)get_block(sizeof(float), n_query * k);

    for (int i = 0; i < n_point * k; i++)
    {
        ds[i] = rand() % 533;
        ds[i] = (float)ds[i] / 10;
    }
    for (int i = 0; i < n_query * k; i++)
    {
        qs[i] = rand() % 800;
        qs[i] = (float)qs[i] / 10;
    }
    sprintf(fname, "%s.ds", filename);
    crea_file(fname, ds, k, n_point);

    sprintf(fname, "%s.qs", filename);
    crea_file(fname, qs, k, n_query);

    //leggi("prova.ds");
    //leggi("prova.qs");

    free(ds);
    free(qs);
}

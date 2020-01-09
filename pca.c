#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

#define MATRIX float *
#define KDTREE float * // modificare con il tipo di dato utilizzato

// void pca(struct params *);
void centraMatrice(MATRIX, int, int);
float calcolaT(float *, int, int, int);
void prodottoMatrice(float *, int, int, MATRIX, int, int, float *, int, int, int, int, int);
void prodottoMatriceTrasp(float *, int, int, MATRIX, int, int, float *, int, int,int, int, int, int);
float norma(float *, int, int, int);
void dividi(float *, int, int, int, float);
void aggiornaDataset(MATRIX, int, int, float *, int, int, float *, int, int, int, int, int);
// void aggiornaDataset(MATRIX, int, int, float *, float *, int, int);

typedef struct
{
    char *filename;     //nome del file, estensione .ds per il data set, estensione .qs per l'eventuale query set
    MATRIX ds;          //data set
    MATRIX qs;          //query set
    int n;              //numero di punti del data set
    int k;              //numero di dimensioni del data/query set
    int nq;             //numero di punti del query set
    int h;              //numero di componenti principali da calcolare 0 se PCA non richiesta
    int kdtree_enabled; //1 per abilitare la costruzione del K-d-Tree, 0 altrimenti
    KDTREE kdtree;      //riferimento al K-d-Tree, NULL se costruzione non richiesta
    float r;            //raggio di query, -1 se range query non richieste
    int silent;         //1 per disabilitare le stampe, 0 altrimenti
    int display;        //1 per stampare i risultati, 0 altrimenti
    MATRIX U;           //matrice U restituita dall'algoritmo PCA
    MATRIX V;           //matrice V restituita dall'algoritmo PCA

    //STRUTTURE OUTPUT MODIFICABILI
    int *QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA; //numero di risposte alle query
} params;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void *get_block(int size, int elements)
{
    return _mm_malloc(elements * size, 16);
}

void free_block(void *p)
{
    _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols)
{
    return (MATRIX)get_block(sizeof(double), rows * cols);
}

void dealloc_matrix(MATRIX mat)
{
    free_block(mat);
}

MATRIX load_data(char *filename, int *n, int *k)
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

    MATRIX data = alloc_matrix(rows, cols);
    status = fread(data, sizeof(float), rows * cols, fp);
    fclose(fp);

    *n = rows;
    *k = cols;

    return data;
}

// PROCEDURE ASSEMBLY
// extern void prova(params* input);

void centraMatrice(MATRIX ds, int n, int k)
{
    int i, j;
    float acc, mean;
    for (j = 0; j < k; j++)
    {
        acc = 0;
        for (i = 0; i < n; i++)
        {
            acc += ds[i * k + j];
        }
        mean = acc / n;
        // printf("\nmean = %f ", mean);
        i = 0;
        for (i = 0; i < n; i++)
        {
            ds[i * k + j] = ds[i * k + j] - mean;
        }
    }
}

float calcolaT(float *a, int n, int k, int cut)
{
    int i;
    float res = 0;
    for (i = 0; i < n; i++)
    {
        res += powf(a[i * k + cut], 2);
    }
    return res;
}

//per dataset transposto rigaA=k, colA=n e rigaB=n;

void prodottoMatriceTrasp(float *result, int rigaRes, int cut, MATRIX ds, int rigaA, int colA, float *vect, int rigaB, int colB, int cut2, int k, int n, int h)
{
    int m, i, j;
    float sum = 0;
    for (m = 0; m < rigaA; m++)
    {
        for (i = 0; i < colB; i++)
        {
            sum = 0;
            for (j = 0; j < rigaB; j++)
            {
                sum += ds[j * k + m] * vect[j * h + cut2];
                printf("\nds = %f , u [cut %d] = %f , sum = %f ", ds[j * k + m],cut2, vect[j * h + cut2], sum);
            }
            result[m * h + cut] = sum;
            printf("V = %f \n", result[m * h + cut]);
        }
    }
    /* Perform multiplication algorithm 
        for( m = 0; m < m1->row; m++ ) {
            for( i = 0; i < m2->column; i++ ) {
                for( j = 0; j < m2->row; j++ ) {

                    (result->matrix)[ m ][ i ] += ( m1->matrix )[ m ][ j ] * ( m2->matrix )[ j ][ i ];

                }
            }
        }
    */
}
void prodottoMatrice(float *result, int rigaRes, int cut, MATRIX ds, int rigaA, int colA, float *vect, int rigaB, int colB, int k, int n, int h)
{
    int m, i, j;
    float sum = 0;
    for (m = 0; m < rigaA; m++)
    {
        for (i = 0; i < colB; i++)
        {
            sum = 0;
            for (j = 0; j < rigaB; j++)
            {
                sum += ds[m * k + j] * vect[j * h + cut];
                // printf("\nds = %f , u = %f , sum = %f ", ds[m * k + j], vect[j * h + cut], sum);
            }
            result[m * h + cut] = sum;
            // printf("V = %f \n", result[m * h + cut]);
        }
    }
}

float norma(float *v, int numRig, int numCol, int cut)
{
    float acc;
    int i;
    for (i = 0; i < numRig; i++)
    {
        acc += pow(v[i * numCol + cut], 2);
    }
    return sqrtf(acc);
}

void dividi(float *v, int numRig, int numCol, int cut, float norm)
{

    int i;
    for (i = 0; i < numRig; i++)
    {
        v[i * numCol + cut] = v[i * numCol + cut] / norm;
    }
}

void aggiornaDataset(MATRIX ds, int numRigDS, int numColDS, float *u, int rigaA, int colA, float *v, int rigaB, int colB, int h, int k, int cut)
{

    // (float *result, int rigaRes, int cut, MATRIX ds, int rigaA, int colA, float *vect, int rigaB, int colB, int k, int n, int h)
    int m, i, j;
    float sum = 0;
    for (m = 0; m < rigaA; m++)
    {
        sum = 0;
        for (i = 0; i < colB; i++)
        {
            // for (j = 0; j < rigaB; j++)
            // {
            sum += u[m * h + cut] * v[i * h + cut];
            // printf("\nu = %f , v= %f , sum = %f ", u[m * h + cut], v[i * h + cut], sum);
            // }
        }
        // printf("  ds = %f \n", ds[m * k + i]);
        ds[m * k + i] -= sum;
        // printf("  ds = %f ", ds[m * k + i]);
    }

    // float res;
    // int m, i, j;
    // for (m = 0; m < numRig; m++)
    // {
    //     for (i = 0; i < numCol; i++)
    //     {
    //         ds[m * numRig + i] -= u[m * h + cut] * v[i * h + cut];

    //     }
    // }
}

/*
*	PCA
* 	=====================
*/
void pca(params *input)
{
    centraMatrice(input->ds, input->n, input->k);
    float theta = 1 * exp(-8);
    // printf("\ntetha= %f ", theta);
    input->V = malloc(input->k * input->h * sizeof(float)); //dimensioni (k x h)
    input->U = malloc(input->n * input->h * sizeof(float)); // dimensioni (n x h)
    if (input->V == NULL || input->U == NULL)
    {
        printf("\nNo MEMORIA");
        exit(1);
    }
    int i;
    for (i = 0; i < input->n; i++)
    {
        input->U[i * input->h] = input->ds[i * input->k];
    }
    int t, t1;
    i = 0;
    float norm, tempV;

    // int n = 6;
    // int k = 4;
    // int size = n * k;
    // // input->k = k;
    // // input->n = n;
    // float *dss = malloc(n * k * sizeof(float));
    // for (int i = 0; i < size; i++)
    //     dss[i] = rand() % 100;
    // printf("\ndataset");
    // for (int i = 0; i < size; i++)
    //     printf("  %f ;", dss[i]);
    // printf("\n");

    // int h = 2;
    // float *v = malloc(k * h * sizeof(float));
    // for (int i = 0; i < k * h; i++)
    //     v[i] = rand() % 100;
    // printf("\nV");
    // for (int i = 0; i < k * h; i++)
    //     printf("  %f ;", v[i]);
    // printf("\n");

    // float *u = malloc(n * h * sizeof(float));
    // for (int i = 0; i < n * h; i++)
    //     u[i] = rand() % 100;
    // printf("\nU");
    // for (int i = 0; i < n * h; i++)
    //     printf("  %f ;", u[i]);
    // printf("\n");

    // aggiornaDataset(dss, n, k, u, n, 1, v, 1, k, h, k, 0);

    // centraMatrice(dss, n, k);

    // printf("\ndataset centrata");
    // for (int i = 0; i < size; i++)
    //     printf("  %f ;", dss[i]);
    // printf("\n");
    // float t = calcolaT(u, n, h, 0);
    // printf("\n t = %f", t);
    // prodottoMatriceTrasp(v, k, 0, dss, k, n, u, n, 1, k, n, h);

    // printf("\nV result prod \n");
    // for (int i = 0; i < k * h; i++)
    //     printf("  %f ;", v[i]);
    // printf("\n");

    // dividi(v, k, h, 0, t);

    // printf("\nV result  div\n");
    // for (int i = 0; i < k * h; i++)
    //     printf("  %f ;", v[i]);
    // printf("\n");

    // float norm = norma(v, k, h, 0);

    // printf("norma = %f ", norm);
    // prodottoMatrice(u, n, 0, dss, n, k, v, k, 1, k, n, h);

    // printf("\nU result \n");
    // for (int i = 0; i < n * h; i++)
    //     printf("  %f ;", u[i]);
    // printf("\n");

    int cont = 0;
    int cut = 0;
    for (i = 0; i < input->h; i++)
    {
        // printf(" inizio iterazione %d ", i);
        do
        {
            
            if(i==0){
                cut = 0;
            }
            else
            {
                cut = i-1;
            }
            

            prodottoMatriceTrasp(input->V, input->k, i, input->ds, input->k, input->n, input->U, input->n, 1, cut, input->k, input->n, input->h);

            t = calcolaT(input->U, input->n, input->h, i);//devo sostituire i con cut
            dividi(input->V, input->k, input->h, i, t);
            norm = norma(input->V, input->k, input->h, i);
            dividi(input->V, input->k, input->h, i, norm);

            prodottoMatrice(input->U, input->n, i, input->ds, input->n, input->k, input->V, input->k, 1, input->k, input->n, input->h);
            tempV = calcolaT(input->V, input->k, input->h, i);
            dividi(input->U, input->n, input->h, i, tempV);
            t1 = calcolaT(input->U, input->n, input->h, i);
            cont++;
            printf("\ncont %d", cont);

        } while (t1 - t >= theta * t1);
        // printf("\n fine iterazione %d ", i);




        aggiornaDataset(input->ds, input->n, input->k, input->U, input->n, 1, input->V, 1, input->k, input->h, input->k, i);
    }
}

int main(int argc, char **argv)
{

    char fname[256];
    int i, j, k;
    clock_t t;
    float time;

    //
    // Imposta i valori di default dei parametri
    //

    params *input = malloc(sizeof(params));

    input->filename = NULL;
    input->h = 0;
    input->kdtree = NULL;
    input->r = -1;
    input->silent = 0;
    input->display = 1;
    input->QA = NULL;
    input->nQA = 0;

    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //

    if (argc <= 1 && !input->silent)
    {
        printf("Usage: %s <data_name> [-pca <h>] [-kdtree [-rq <r>]]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\t-d: display query results\n");
        printf("\t-s: silent\n");
        printf("\t-pca <h>: h-component PCA enabled\n");
        printf("\t-kdtree: kdtree building enabled\n");
        printf("\t-rq <r>: range query search with radius r enabled\n");
        printf("\n");
        exit(0);
    }

    //
    // Legge i valori dei parametri da riga comandi
    //

    int par = 1;
    while (par < argc)
    {
        if (par == 1)
        {
            input->filename = argv[par];
            par++;
        }
        else if (strcmp(argv[par], "-s") == 0)
        {
            input->silent = 1;
            par++;
        }
        else if (strcmp(argv[par], "-d") == 0)
        {
            input->display = 1;
            par++;
        }
        else if (strcmp(argv[par], "-pca") == 0)
        {
            par++;
            if (par >= argc)
            {
                printf("Missing h value!\n");
                exit(1);
            }
            input->h = atoi(argv[par]);
            par++;
        }
        else if (strcmp(argv[par], "-kdtree") == 0)
        {
            input->kdtree_enabled = 1;
            par++;
            if (par < argc && strcmp(argv[par], "-rq") == 0)
            {
                par++;
                if (par >= argc)
                {
                    printf("Missing radius value!\n");
                    exit(1);
                }
                input->r = atof(argv[par]);
                if (input->r < 0)
                {
                    printf("Range query radius must be non-negative!\n");
                    exit(1);
                }
                par++;
            }
        }
        else
        {
            printf("WARNING: unrecognized parameter '%s'!\n", argv[par]);
            par++;
        }
    }

    //
    // Legge i dati e verifica la correttezza dei parametri
    //

    if (input->filename == NULL || strlen(input->filename) == 0)
    {
        printf("Missing input file name!\n");
        exit(1);
    }

    sprintf(fname, "%s.ds", input->filename);
    input->ds = load_data(fname, &input->n, &input->k);

    if (input->h < 0)
    {
        printf("Invalid value of PCA parameter h!\n");
        exit(1);
    }
    if (input->h > input->k)
    {
        printf("Value of PCA parameter h exceeds data set dimensions!\n");
        exit(1);
    }

    if (input->r >= 0)
    {
        sprintf(fname, "%s.qs", input->filename);
        int k;
        input->qs = load_data(fname, &input->nq, &k);
        if (input->k != k)
        {
            printf("Data set dimensions and query set dimensions are not compatible!\n");
            exit(1);
        }
    }

    //
    // Visualizza il valore dei parametri
    //

    if (!input->silent)
    {
        printf("Input file name: '%s'\n", input->filename);
        printf("Data set size [n]: %d\n", input->n);
        printf("Number of dimensions [k]: %d\n", input->k);
        if (input->h > 0)
        {
            printf("PCA search enabled\n");
            printf("Number of principal components [h]: %i\n", input->h);
        }
        else
        {
            printf("PCA search disabled\n");
        }
        if (input->kdtree_enabled)
        {
            printf("Kdtree building enabled\n");
            if (input->r >= 0)
            {
                printf("Range query search enabled\n");
                printf("Range query search radius [r]: %f\n", input->r);
            }
            else
            {
                printf("Range query search disabled\n");
            }
        }
        else
        {
            printf("Kdtree building disabled\n");
        }
    }

    //
    // Calcolo PCA
    //

    if (input->h > 0)
    {
        printf("\nINIZIO PCA\n");
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t) / CLOCKS_PER_SEC;
        // sprintf(fname, "%s.U", input->filename);
        // save_data(fname, input->U, input->n, input->h);
        // sprintf(fname, "%s.V", input->filename);
        // save_data(fname, input->V, input->k, input->h);
    }
    else
        time = -1;

    if (!input->silent)
        printf("\nPCA time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    return 0;
}
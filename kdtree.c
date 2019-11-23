#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

#define MATRIX float *
#define KDTREE float * // modificare con il tipo di dato utilizzato

int partition2(MATRIX, int *, int, int, int, int);
void MyprintArray(MATRIX, int *, int, int, int);
int partition(MATRIX, int *, int, int, int, int);
void swap(int *, int *);
void quicksort(int *, MATRIX, int, int, int, int, int);
void buildTree(MATRIX, int, int, int);

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

typedef struct
{
    float *point;
    float h_min, h_max;
    struct kdtree_node *left, *right;
} kdtree_node;

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

    printf("\ncolonne= %d\trighe= %d\n", cols, rows);

    MATRIX data = alloc_matrix(rows, cols);
    status = fread(data, sizeof(float), rows * cols, fp);
    fclose(fp);

    *n = rows;
    *k = cols;

    return data;
}

void MyprintArray(MATRIX ds, int *a, int size, int k, int cut)
{
    printf("\nInizio a stampare \n");
    for (int i = 0; i < size; i++)
    {
        printf("i: %d - %f, ", a[i], ds[a[i] * k + cut]);
        if (i % 10 == 0 && i > 0)
            printf("\n");
    }
    printf("\n");
}

void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(MATRIX dataset, int *indexSorted, int low, int high, int k, int cut)
{

    float pivot = dataset[(high - 1) * k + cut];
    int i = (low - 1);

    for (int j = low; j < high; j++)
    {
        if (dataset[j * k + cut] < pivot)
        {
            i++;
            //swap
            indexSorted[i] = i;
            indexSorted[j] = j;
        }
        else
        {
            indexSorted[j] = j;
        }
    }
    //final swap
    indexSorted[high] = i + 1;
    indexSorted[i + 1] = high;
    return (i + 1);
}

int partition2(MATRIX dataset, int *indexSorted, int low, int high, int k, int cut)
{

    float pivot = dataset[indexSorted[high] * k + cut];
    int i = (low - 1);

    for (int j = low; j < high; j++)
    {
        if (dataset[indexSorted[j] * k + cut] < pivot)
        {
            i++;
            swap(&indexSorted[i], &indexSorted[j]);
        }
    }
    swap(&indexSorted[i + 1], &indexSorted[high]);
    return (i + 1);
}

void quicksort(int *indexSorted, MATRIX dataset, int cut, int k, int low, int high, int flagg)
{
    if (low < high)
    {

        int pi = -1;
        if (flagg == 0)
        {
            //solo la prima volta che chiamo il metodo devo utilizzare il dataset
            //riempio array indice con gli indici ordinati partialmente (solo 1 itarazione)
            pi = partition(dataset, indexSorted, low, high, k, cut);
        }
        else
        {
            //fatta la prima chiamata ricorsiva
            //uso come indici quelli memorizzati nell'array parzialemente ordinato
            pi = partition2(dataset, indexSorted, low, high, k, cut);
        }
        quicksort(indexSorted, dataset, cut, k, low, pi - 1, 1);
        quicksort(indexSorted, dataset, cut, k, pi + 1, high, 1);
    }
}

// struct kdtree_node* buildTree(MATRIX ds, int liv, int size, int k)
void buildTree(MATRIX ds, int liv, int size, int k)
{
    if (size == 0)
    {
        printf("ERRORE SIZE NULLA");
        return;
    }

    int cut = liv % k;                             //variabile di cut
    int *indexSorted = malloc(size * sizeof(int)); //vettore che conterra indice riga dei punti ordinati
    quicksort(indexSorted, ds, cut, k, 0, size, 0);

    //serve per stampare l'array di indici e i punti associati
    MyprintArray(ds, indexSorted, size, k, cut);
    printf("\nval min= %f\tval max= %f\t valore Mediano= %f\n", ds[indexSorted[0] * k + cut], ds[indexSorted[size - 1] * k + cut], ds[indexSorted[size / 2] * k + cut]);

    //PSEUDO CODE
    //ordinare colonna dataset[ki] con un metodo  sort
    //ordinare gli indici
    //prendere punto mediano P
    //prendere indici datasset <P e >= P
    // kdtree_node *node = malloc(sizeof(kdtree_node));
    // node->point = P mediano
    // node->h_min = bound min
    // node->h_max = bound max
    // node->left = buildTree(ds, liv + 1, size - 1, k);
    // node->right = buildTree(ds, liv + 1, size - 1, k);

    // return node;
}

/*
*	K-d-Tree
* 	======================
*/
void kdtree(params *input)
{
    printf("INizio kdtree");
    printf("\ndataset size%d, dataset k%d\n", input->n, input->k);

    buildTree(input->ds, 0, input->n, input->k);
}

int main(int argc, char const *argv[])
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

    t = clock();
    kdtree(input);
    t = clock() - t;
    time = ((float)t) / CLOCKS_PER_SEC;
    // if (input->kdtree)
    // {
    //     t = clock();
    //     kdtree(input);
    //     t = clock() - t;
    //     time = ((float)t) / CLOCKS_PER_SEC;
    // }
    // else
    //     time = -1;

    return 0;
}

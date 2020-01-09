#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <malloc.h>

#define MATRIX float *
#define KDTREE struct kdtree_node * // modificare con il tipo di dato utilizzato

void swap(int *, int *);
int partition(MATRIX, int *, int, int, int, int);
void MyprintArray(MATRIX, int *, int, int, int, int);
void quicksort(MATRIX, int *, int, int, int, int);
int findMedian(MATRIX, int *, int, int, int, int);
struct kdtree_node *buildTreeRoot(MATRIX, int *, int, int, int);
struct kdtree_node *buildTree(MATRIX, int *, int, int, int, int, int, int, float *);
float *findRegion(MATRIX, int, int);
float euclidean_distance(MATRIX qs, int id_qs, MATRIX ds, int id_ds, int k);
int *rangequery(MATRIX ds, struct kdtree_node *tree, MATRIX qs, int id_qs, float r, int k, int n, int *list, int *nQA);

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
    MATRIX region;      //regione indicizzata da root formata da k coppie (h_min, h_max)

    //STRUTTURE OUTPUT MODIFICABILI
    int *QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA; //numero di risposte alle query
} params;

struct kdtree_node
{
    float medianCoordinate; //coordinata cut per il punto mediano
    int indexMedianPoint;
    int numPoint;
    float *region;
    struct kdtree_node *left, *right;
};

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

/*
* 
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
* 
*/
void save_data(char *filename, void *X, int n, int k)
{
    FILE *fp;
    int i;
    fp = fopen(filename, "wb");
    if (X != NULL)
    {
        fwrite(&n, 4, 1, fp);
        fwrite(&k, 4, 1, fp);
        for (i = 0; i < n; i++)
        {
            fwrite(X, 4, k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += 4 * k;
        }
    }
    fclose(fp);
}

/*
*	PCA
* 	=====================
*/
void pca(params *input)
{

    // -------------------------------------------------
    // Codificare qui l'algoritmo PCA
    // -------------------------------------------------
    // prova(input);
    // Calcola le matrici U e V
    // -------------------------------------------------
}

void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(MATRIX dataset, int *indexSorted, int low, int high, int k, int cut)
{
    float pivot = dataset[indexSorted[high - 1] * k + cut];
    int i = (low - 1);
    int j;
    for (j = low; j < high - 1; j++)
    {
        if (dataset[indexSorted[j] * k + cut] < pivot)
        {
            i++;
            swap(&indexSorted[i], &indexSorted[j]);
        }
    }

    swap(&indexSorted[i + 1], &indexSorted[high - 1]);
    return (i + 1);
}

void quicksort(MATRIX dataset, int *indexSorted, int cut, int k, int low, int high)
{
    if (low < high)
    {
        int indexPivot = -1;

        indexPivot = partition(dataset, indexSorted, low, high, k, cut);

        quicksort(dataset, indexSorted, cut, k, low, indexPivot);
        quicksort(dataset, indexSorted, cut, k, indexPivot + 1, high);
    }
}

int findMedian(MATRIX dataset, int *indexSorted, int start, int indexMedian, int k, int cut)
{
    float median = dataset[indexSorted[indexMedian] * k + cut];
    int i;
    for (i = indexMedian - 1; i >= start; i--)
    {
        if (dataset[indexSorted[i] * k + cut] < median)
        {
            return i + 1;
        }
    }
    return indexMedian;
}

struct kdtree_node *buildTree(MATRIX ds, int *indexSorted, int liv, int start, int end, int numEle, int k, int type, float *regionp)
{
    if (numEle == 0)
        return NULL;

    int cut = liv % k; //variabile di cut per indice colonna da usare
    struct kdtree_node *node = (struct kdtree_node *)malloc(sizeof(struct kdtree_node));
    node->region = malloc(2 * k * sizeof(float));
    for (int i = 0; i < (2 * k); i++)
        node->region[i] = regionp[i];
    node->region[2 * (cut - 1)] = ds[indexSorted[start] * k + (cut - 1)];
    node->region[2 * (cut - 1) + 1] = ds[indexSorted[end - 1] * k + (cut - 1)];

    node->numPoint = numEle;
    quicksort(ds, indexSorted, cut, k, start, end);

    int indexMedian = findMedian(ds, indexSorted, start, start + ((end - 1 - start) / 2), k, cut);

    node->medianCoordinate = ds[indexSorted[indexMedian] * k + cut]; //valore di coordinata del punto mediano
    node->indexMedianPoint = indexSorted[indexMedian];               //indice del punto mediano nel dataset

    int numEleSx = indexMedian - start;
    int numEleDx = end - indexMedian - 1;

    if (numEleSx == 0 && numEleDx == 0)
    {
        node->left = NULL;
        node->right = NULL;
        return node;
    }
    else if (numEleSx == 0)
    {
        node->left = NULL;
        node->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1, node->region);
        return node;
    }
    else if (numEleDx == 0)
    {
        node->right = NULL;
        node->left = buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0, node->region);
        return node;
    }
    else
    {
        node->left = buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0, node->region);
        node->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1, node->region);
        return node;
    }
}

struct kdtree_node *buildTreeRoot(MATRIX ds, int *indexSorted, int liv, int end, int k)
{
    if (end <= 0)
    {
        printf("\nDATASET NULLO\n");
        return NULL;
    }
    int cut = liv % k; //variabile di cut per indice colonna da usare

    quicksort(ds, indexSorted, cut, k, 0, end);

    int indexMedian = findMedian(ds, indexSorted, 0, (end - 1) / 2, k, cut);

    struct kdtree_node *root = (struct kdtree_node *)malloc(sizeof(struct kdtree_node));
    root->medianCoordinate = ds[indexSorted[indexMedian] * k + cut]; //valore di coordinata del punto mediano
    root->indexMedianPoint = indexSorted[indexMedian];               //indice del punto mediano nel dataset
    root->region = findRegion(ds, end, k);
    root->numPoint = end;
    int numEleSx = indexMedian;
    int numEleDx = end - indexMedian - 1;

    root->left = buildTree(ds, indexSorted, liv + 1, 0, indexMedian, numEleSx, k, 0, root->region);
    root->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1, root->region);

    return root;
}

float *findRegion(MATRIX ds, int n, int k)
{
    float *region = (float *)malloc(k * 2 * sizeof(float));
    float h_min, h_max;
    int j, i;
    for (j = 0; j < k; j++)
    {
        h_min = ds[j];
        h_max = ds[j];
        for (i = 0; i < n; i++)
        {
            if (h_max < ds[i * k + j])
                h_max = ds[i * k + j];
            if (h_min > ds[i * k + j])
                h_min = ds[i * k + j];
        }
        region[2 * j] = h_min;
        region[(2 * j) + 1] = h_max;
    }
    return region;
}

int indexList = 0;

void kdtree(params *input)
{

    printf("\nInizio kdtree");
    int *indexSorted = (int *)malloc(input->n * sizeof(int)); //vettore che conterra indice riga dei punti ordinati

    if (indexSorted == NULL)
    {
        printf("\nNO MEMORIA\n");
        exit(1);
    }
    int i;
    for (i = 0; i < input->n; i++)
    {
        indexSorted[i] = i;
    }
    // si può usare memset
    // memset(indexSorted,0,input->n*sizeof(int));

    input->kdtree = buildTreeRoot(input->ds, indexSorted, 0, input->n, input->k);

    // // printTree(input->kdtree);

    //bisogna liberare la memoria
    free(indexSorted);
    printf("\nfine KDTREE");

    // free(region);
}

float distance(float *h, MATRIX qs, int id_qs, int k)
{
    int i = 0;
    float *p = malloc(k * sizeof(float)); //serve per creare il punto sul bordo
    for (; i < k; i++)
    {
        if (qs[k * id_qs + i] <= h[i * 2])
            p[i] = h[i * 2];
        else if (qs[k * id_qs + i] >= h[(i * 2) + 1])
            p[i] = h[(i * 2) + 1];
        else
            p[i] = qs[k * id_qs + i];
    } //for
    float temp = euclidean_distance(qs, id_qs, p, 0, k);
    free(p);
    return temp;
}

float euclidean_distance(MATRIX qs, int id_qs, MATRIX ds, int id_ds, int k)
{
    float somma = 0;

    int i = 0;
    for (; i < k; i++)
    {
        somma = somma + (qs[id_qs * k + i] - ds[id_ds * k + i]) * (qs[id_qs * k + i] - ds[id_ds * k + i]);
    }
    return sqrtf(somma);
}

int *rangequery(MATRIX ds, struct kdtree_node *tree, MATRIX qs, int id_qs, float r, int k, int n, int *list, int *nQA)
{
    if (tree == NULL || distance(tree->region, qs, id_qs, k) > r){
        // printf("stop  ");
        return list;
    }

    if (euclidean_distance(qs, id_qs, ds, tree->indexMedianPoint, k) <= r)
    {
        list[id_qs * n + indexList] = tree->indexMedianPoint;
        indexList++;
        *nQA = *nQA + 1;
    };
    if (tree->left != NULL)
    {
        rangequery(ds, tree->left, qs, id_qs, r, k, n, list, nQA);
        // return list;
    }
    if (tree->right != NULL)
    {
        rangequery(ds, tree->right, qs, id_qs, r, k, n, list, nQA);
        // return list;
    }

    return list;
} //rangequery

void range_query(params *input)
{
    printf("\n inizio QUERY");

    input->QA = (int *)malloc(input->n * input->nq * sizeof(int));
    if (input->QA == NULL)
    {
        printf("\nNO MEMORIA");
        exit(1);
    }
    int i;
    for (i = 0; i < input->nq; i++)
    {

        rangequery(input->ds, input->kdtree, input->qs, i, input->r, input->k, input->n, input->QA, &input->nQA);
        input->QA[i * input->n + indexList] = -1;
        indexList = 0;
    }
}

void stampa(float *prova, int size)
{
    for (int i = 0; i < size - 1; i += 2)
    {
        printf("%f,%f\n", prova[i], prova[i + 1]);
    }
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
    input->r = 0;
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

    // t = clock();
    // kdtree(input);
    // t = clock() - t;
    // time = ((float)t) / CLOCKS_PER_SEC;
    // printf("\n\ntime= %f seconds\n", time);

    // t = clock();
    // range_query(input);
    // t = clock() - t;
    // time = ((float)t) / CLOCKS_PER_SEC;
    // printf("\n\ntime= %f seconds\n", time);

    if (input->h > 0)
    {
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t) / CLOCKS_PER_SEC;
        sprintf(fname, "%s.U", input->filename);
        save_data(fname, input->U, input->n, input->h);
        sprintf(fname, "%s.V", input->filename);
        save_data(fname, input->V, input->k, input->h);
    }
    else
        time = -1;

    if (!input->silent)
        printf("\nPCA time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Costruzione K-d-Tree
    //

    // if (input->kdtree)
    // {
    t = clock();
    kdtree(input);
    t = clock() - t;
    time = ((float)t) / CLOCKS_PER_SEC;
    // }
    // else
    //     time = -1;
    if (!input->silent)
        printf("\nIndexing time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Range query search
    //

    if (input->r >= 0)
    {
        t = clock();
        range_query(input);
        t = clock() - t;
        time = ((float)t) / CLOCKS_PER_SEC;
    }
    else
        time = -1;
    if (!input->silent)
        printf("\nQuerying time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Salva il risultato delle query
    // da modificare se si modifica il formato delle matrici di output
    //

    if (input->r >= 0)
    {
        if (!input->silent && input->display)
        {
            //NB: il codice non assume che QA sia ordinata per query, in caso lo sia ottimizzare il codice
            printf("\nQuery Answer:\n");
            int flag, j;
            for (i = 0; i < input->nq; i++)
            {
                for (j = 0; j < input->n; j++)
                {
                    if (input->QA[i * input->n + j] == -1)
                    {
                        // printf("]");
                        break;
                    }
                    if (flag == 0)
                    {
                        printf("\nid_q %d: [", i);
                        flag = 1;
                    }
                    printf("%d, ", input->QA[i * input->n + j]);
                }
                flag = 0;
            }
            printf("\n");
        }
        sprintf(fname, "%s.qa", input->filename);
        // save_data(fname, input->QA, input->nQA, 2);
    }

    if (!input->silent)
        printf("\nDone.\n");

    return 0;
}
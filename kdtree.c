#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <malloc.h>

#define MATRIX float *
#define KDTREE float * // modificare con il tipo di dato utilizzato

void swap(int *, int *);
int partition(MATRIX, int *, int, int, int, int);
int partition1(MATRIX, int *, int *, int, int, int, int);
int partition2(MATRIX, int *, int, int, int, int);
void MyprintArray(MATRIX, int *, int, int, int);
void quicksort(MATRIX, int *, int, int, int, int, int);
void quicksort2(MATRIX, int *, int *, int, int, int, int, int);
int findMedian(MATRIX, int *, int, int, int);
// struct kdtree_node *buildTreeRoot(MATRIX, int *, int, int, int);
// struct kdtree_node *buildTree(MATRIX, int *, int *, int, int, int);
void buildTreeRoot(MATRIX, int *, int, int, int);
void buildTree(MATRIX, int *, int *, int, int, int);
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
        // if (i == 100)
        //     break;
    }
    printf("\n");
}

void printArray(int *a, int size)
{
    printf("\nInizio a stampare \n");
    for (int i = 0; i < size; i++)
    {
        printf("i: %d - %d ", i, a[i]);
        if (i % 10 == 0 && i > 0)
            printf("\n");
        if (i == 100)
            break;
    }
    printf("\n");
}

void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}
/* 
*   partition viene usato solo per il nodo root alla prima chiamata
*   per riempire il vettore degli indici con un ordinamento parziale
*/
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
    indexSorted[high - 1] = i + 1;
    indexSorted[i + 1] = high - 1;
    return (i + 1);
}

/*
*   partition1 viene usato alla prima iterazione per tutti i nodi del kd tree che non sono root
*   poichè consideriamo solo la parte di punti minori o maggiori della mediana in cui memorizziamo
*   l'indice in indexSorted, mentre in newIndexSorted ci saranno gli indici dei punti ordinati secondo
*   un'altra coordinata k
*/
int partition1(MATRIX dataset, int *indexSorted, int *newIndexSorted, int low, int high, int k, int cut)
{

    float pivot = dataset[indexSorted[(high - 1)] * k + cut];
    int i = (low - 1);

    for (int j = low; j < high; j++)
    {
        if (dataset[indexSorted[j] * k + cut] < pivot)
        {
            i++;
            //swap
            newIndexSorted[i] = i;
            newIndexSorted[j] = j;
        }
        else
        {
            newIndexSorted[j] = j;
        }
    }
    //final swap
    newIndexSorted[high - 1] = i + 1;
    newIndexSorted[i + 1] = high - 1;
    return (i + 1);
}

/*  
*   partition2 viene usato dalla seconda iterazione in cui si prende l'indice dei punti da ordinare
*   nel vettore indexSorted
*/
int partition2(MATRIX dataset, int *indexSorted, int low, int high, int k, int cut)
{
    // printf("\npartition2   low= %d, high= %d  cut= %d", low, high, cut);
    // printf("\n indexSorted[hight-1]= %d, high= %d", indexSorted[high - 1], high);
    // printf("\n indexSorted[hight-1]= %d, high= %d", indexSorted[high], high);

    float pivot = dataset[indexSorted[high - 1] * k + cut];
    // printf("\narr val pivot= %f, indexSorted[hight-1]= %d, high= %d", pivot, indexSorted[high - 1], high);
    // printf("\n indexSorted[hight-1]= %d, high= %d", indexSorted[high - 1], high);
    // printf("\n indexSorted[hight-1]= %d, high= %d\n", indexSorted[high], high);
    int i = (low - 1);
    for (int j = low; j < high; j++)
    {
        if (dataset[indexSorted[j] * k + cut] < pivot)
        {
            i++;
            // printf("\tindex= %d val= %f", indexSorted[j], dataset[indexSorted[j] * k + cut]);
            swap(&indexSorted[i], &indexSorted[j]);
        }
    }
    // printf("\n indexSorted= %p", indexSorted);
    // printf("\n\n");
    // printf("\ni= %d", i);

    swap(&indexSorted[i + 1], &indexSorted[high - 1]);
    return (i + 1);
}

void quicksort2(MATRIX dataset, int *indexSorted, int *newIndexSorted, int cut, int k, int low, int high, int flagg)
{
    if (low < high)
    {
        // printf("\nflag %d, low= %d, high= %d", flagg, low, high);
        int indexPivot = -1;
        if (flagg == 0)
        {
            // printf("\nprima del sort 1 parzione\n");
            // MyprintArray(dataset, indexSorted, high, k, cut);
            //solo per il nodo root e una sola chiamata
            //riempio array indice con gli indici ordinati partialmente (solo 1 iterazione)
            indexPivot = partition1(dataset, indexSorted, newIndexSorted, low, high, k, cut);
            // printf("\n1 partizione figlio e index pivot= %d", indexPivot);
            // MyprintArray(dataset, newIndexSorted, high, k, cut);
            // printf("\nindexSorted:");
            // printArray(indexSorted, high);
            // printf("\nnew Index array");
            // printArray(newIndexSorted, high);
        }
        else
        {
            // printf("\nseconda partizione figlio");
            // printf("\n indexSorted= %p", indexSorted);
            // printf("\n newIdexSorted= %p", newIndexSorted);
            //  printf("\nindexSorted:");
            // printArray(indexSorted, high);
            // printf("\nnew Index array");
            // printArray(newIndexSorted, high);
            // printf("\n\n");

            //fatta la prima chiamata ricorsiva
            //uso come indici quelli memorizzati nell'array parzialemente ordinato
            indexPivot = partition2(dataset, newIndexSorted, low, high, k, cut);
        }

        quicksort2(dataset, indexSorted, newIndexSorted, cut, k, low, indexPivot, 1);
        quicksort2(dataset, indexSorted, newIndexSorted, cut, k, indexPivot + 1, high, 1);
    }
}

void quicksort(MATRIX dataset, int *indexSorted, int cut, int k, int low, int high, int flagg)
{
    if (low < high)
    {

        int indexPivot = -1;
        if (flagg == 0)
        {
            //solo per il nodo root e una sola chiamata
            //riempio array indice con gli indici ordinati partialmente (solo 1 iterazione)
            indexPivot = partition(dataset, indexSorted, low, high, k, cut);
            flagg = 1;
        }
        else
        {
            //fatta la prima chiamata ricorsiva
            //uso come indici quelli memorizzati nell'array parzialemente ordinato
            indexPivot = partition2(dataset, indexSorted, low, high, k, cut);
        }

        quicksort(dataset, indexSorted, cut, k, low, indexPivot - 1, 1);
        quicksort(dataset, indexSorted, cut, k, indexPivot + 1, high, 1);
    }
}

int findMedian(MATRIX dataset, int *indexSorted, int indexMedian, int k, int cut)
{
    float median = dataset[indexSorted[indexMedian] * k + cut];
    for (int i = indexMedian - 1; i >= 0; i--)
    {
        if (dataset[indexSorted[i] * k + cut] < median)
        {
            printf("\nfind median  val median= %f indexMedianVecchio= %d", median, indexMedian);
            printf("\nNUOVO  val median= %f  indexMedian= %d valMedian[index-1]= %f  valMedian[index+1]= %f", median, i + 1, dataset[indexSorted[i] * k + cut], dataset[indexSorted[i + 2] * k + cut]);
            return i + 1;
        }
    }
    return -1;
}

/*
*   buildTree server per costruire tutti i nodi del kdtree
*/
// struct kdtree_node *buildTree(MATRIX ds, int *indexSorted, int *newIndexSorted, int liv, int size, int k)
void buildTree(MATRIX ds, int *indexSorted, int *newIndexSorted, int liv, int size, int k)
{
    if (size <= 0)
    {
        printf("\nSIZE NULLA\n");
        // return;
    }
    printf("\nCOSTRIUSCO FIGLIO liv= %d, size= %d, \n", liv, size);

    int cut = liv % k; //variabile di cut per indice colonna da usare

    if (cut == 2)
    {
        printf("\nFINE\n");
        return;
    }
    MyprintArray(ds, indexSorted, size, k, cut - 1);
    quicksort2(ds, indexSorted, newIndexSorted, cut, k, 0, size, 0);

    // free(indexSorted);
    // if (indexSorted == NULL)
    //     printf("\ndeallocato indexSOrted");
    //POSSO DEALLOCARE QUI INDEXSORTED PERCHE NON LO UTILIZZO PIU
    int indexMedian = findMedian(ds, newIndexSorted, size / 2, k, cut);
    // int indexMedian = newIndexSorted[size / 2];
    //serve per stampare l'array di indici e i punti associati

    // MyprintArray(ds, newIndexSorted, size, k, cut - 1);
    printf("\nval min= %f\tval max= %f\t valore Mediano= %f", ds[newIndexSorted[0] * k + cut], ds[newIndexSorted[size - 1] * k + cut], ds[newIndexSorted[indexMedian] * k + cut]);
    kdtree_node *node = malloc(sizeof(kdtree_node));
    //SUGGERIMENTO: potrebbe bastare memorizzare solo una coordinata invece di tutto il punto con le k coordinate???
    node->point = &ds[indexMedian];             //in point ci sarà l'indirizzo della prima cella del mediano
    node->h_min = ds[newIndexSorted[0]];        //valore di coordinata più piccola
    node->h_max = ds[newIndexSorted[size - 1]]; //valore di coordinata più grande

    int newSizeSx = indexMedian + 1;
    int newSizeDx = size - newSizeSx;
    int *newIndexSortedSx = malloc(newSizeSx * sizeof(int));
    int *newIndexSortedDx = malloc(newSizeDx * sizeof(int));
    printf("\nsize iniziale= %d sizeSX= %d  sizeDX= %d", size, newSizeSx, newSizeDx);
    if (newIndexSortedDx == NULL | newIndexSortedSx == NULL)
    {
        printf("\nNO MEMORIA FIGLI\n");
        exit(1);
    }
    //indexSorted = riferimento prima cella della prima metà di punti
    buildTree(ds, newIndexSorted, newIndexSortedSx, liv + 1, newSizeSx, k);
    buildTree(ds, newIndexSorted + (newSizeSx), newIndexSortedDx, liv + 1, newSizeDx, k);
    // node->left = buildTree(ds, newIndexSorted, newIndexSortedSx, liv + 1, (size / 2), k);
    //indexSorted + ((size / 2) + 1) = riferimento prima cella della seconda metà di punti
    // node->right = buildTree(ds, newIndexSorted + ((size / 2) + 1), newIndexSortedSx, liv + 1, (size / 2), k);

    // return node;
}

/*
*   buildTreeRoot server per costruire solo la radice del hdtree
*/
// struct kdtree_node *buildTreeRoot(MATRIX ds, int *indexSorted, int liv, int size, int k)
void buildTreeRoot(MATRIX ds, int *indexSorted, int liv, int size, int k)
{
    if (size <= 0)
    {
        printf("\nSIZE NULLA\n");
        // return;
    }
    printf("\nCOSTRIUSCO RADICE\n");
    int cut = liv % k; //variabile di cut per indice colonna da usare
    quicksort(ds, indexSorted, cut, k, 0, size, 0);

    int indexMedian = findMedian(ds, indexSorted, size / 2, k, cut);
    // int indexMedian = indexSorted[size / 2];

    //serve per stampare l'array di indici e i punti associati
    // MyprintArray(ds, indexSorted, size, k, cut);
    printf("\nval min= %f\tval max= %f\t valore Mediano= %f\n", ds[indexSorted[0] * k + cut], ds[indexSorted[size - 1] * k + cut], ds[indexSorted[indexMedian] * k + cut]);
    kdtree_node *root = malloc(sizeof(kdtree_node));
    //SUGGERIMENTO: potrebbe bastare memorizzare solo una coordinata invece di tutto il punto con le k coordinate???
    root->point = &ds[indexMedian];          //in point ci sarà l'indirizzo della prima cella del mediano
    root->h_min = ds[indexSorted[0]];        //valore di coordinata più piccola
    root->h_max = ds[indexSorted[size - 1]]; //valore di coordinata più grande

    //PSEUDO CODE
    //ordinare colonna dataset[ki] con un metodo  sort
    //ordinare gli indici
    //prendere punto mediano P
    //prendere indici datasset <P e >= P
    // kdtree_node *node = malloc(sizeof(kdtree_node));
    // node->point = P mediano
    // node->h_min = bound min
    // node->h_max = bound max
    int newSizeSx = indexMedian + 1;
    int newSizeDx = size - newSizeSx;
    int *newIndexSortedSx = malloc(newSizeSx * sizeof(int));
    int *newIndexSortedDx = malloc(newSizeDx * sizeof(int));
    printf("\nsize iniziale= %d sizeSX= %d  sizeDX= %d", size, newSizeSx, newSizeDx);
    if (newIndexSortedDx == NULL | newIndexSortedSx == NULL)
    {
        printf("\nNO MEMORIA FIGLI\n");
        exit(1);
    }
    buildTree(ds, indexSorted, newIndexSortedSx, liv + 1, newSizeSx, k);
    buildTree(ds, indexSorted + (newSizeSx), newIndexSortedDx, liv + 1, newSizeDx, k);
    // buildTree(ds, indexSorted, newIndexSortedSx, liv + 1, (size / 2), k);
    // buildTree(ds, indexSorted + (size / 2), newIndexSortedDx, liv + 1, (size / 2), k);
    //indexSorted = riferimento prima cella della prima metà di punti
    // root->left = buildTree(ds, indexSorted, newIndexSortedSx, liv + 1, (size / 2) - 1, k);
    //indexSorted + ((size / 2) + 1) = riferimento prima cella della seconda metà di punti
    // root->right = buildTree(ds, indexSorted + ((size / 2) + 1), newIndexSortedSx, liv + 1, (size / 2) - 1, k);

    // return root;
}

/*
*	K-d-Tree
* 	======================
*/
void kdtree(params *input)
{
    printf("INizio kdtree");
    printf("\ndataset size%d, dataset k%d\n", input->n, input->k);
    int *indexSorted = (int *)malloc(input->n * sizeof(int)); //vettore che conterra indice riga dei punti ordinati
    if (indexSorted == NULL)
    {
        printf("\nNO MEMORIA\n");
        exit(1);
    }
    buildTreeRoot(input->ds, indexSorted, 0, input->n, input->k);
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

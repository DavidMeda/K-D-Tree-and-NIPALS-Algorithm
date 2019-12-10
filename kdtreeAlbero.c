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
struct kdtree_node *buildTree(MATRIX, int *, int, int, int, int, int, int);

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

struct kdtree_node
{
    float median; //median coordinata
    float h_min, h_max;
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

void MyprintArray(MATRIX ds, int *a, int start, int size, int k, int cut)
{
    printf("\nInizio a stampare \n");
    for (int i = start; i < size; i++)
    {
        printf("i:%d - %f, ", a[i], ds[a[i] * k + cut]);
        if (i % 10 == 0 && i > 0)
            printf("\n");
        if (i == 100)
            break;
    }
    printf("\n");
}

void printArray(int *a, int size)
{
    printf("\nInizio a stampare \n");
    for (int i = 0; i < size; i++)
    {
        printf("i: %d - %d, ", i, a[i]);
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

int partition(MATRIX dataset, int *indexSorted, int low, int high, int k, int cut)
{
    float pivot = dataset[indexSorted[high - 1] * k + cut];
    int i = (low - 1);
    for (int j = low; j < high - 1; j++)
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
    for (int i = indexMedian - 1; i >= start; i--)
    {
        if (dataset[indexSorted[i] * k + cut] < median)
        {
            return i + 1;
        }
    }
    return indexMedian;
}
//contatore per il debug
int cont = 1;

/*
*   buildTree serve per costruire tutti i nodi del kdtree
*   il valore end deve essere escluso indici vanno da [start, end) quindi end è escluso l'ultimo elemento si trova a end-1
*/
struct kdtree_node *buildTree(MATRIX ds, int *indexSorted, int liv, int start, int end, int numEle, int k, int type)
{

    int cut = liv % k; //variabile di cut per indice colonna da usare
    cont++;
    // type = 0 nodo sinistro, type = 1 node destro (serve per il debug)
    // if (type == 0)
    //     printf("\n\nCOSTRIUSCO FIGLIO SINISTRO liv= %d, start= %d end= %d, cut= %d, numEle= %d, cont= %d", liv, start, end, cut, numEle, cont);
    // else
    //     printf("\n\nCOSTRIUSCO FIGLIO DESTRO liv= %d, start= %d end= %d, cut= %d, numEle= %d, cont= %d", liv, start, end, cut, numEle, cont);

    //serve per il debug
    // if (cut == 5)
    // {
    //     // printf("\nFINE\n");
    //     return NULL;
    // }
    quicksort(ds, indexSorted, cut, k, start, end);

    int indexMedian = findMedian(ds, indexSorted, start, start + ((end - 1 - start) / 2), k, cut);

    // printf("\nvalore MEDIANO= %f\t indexMedian=%d\t val puntoMin= %f\t val puntMax= %f ", ds[indexSorted[indexMedian] * k + cut], indexMedian, ds[indexSorted[start] * k + cut], ds[indexSorted[end - 1] * k + cut]);

    //SUGGERIMENTO: potrebbe bastare memorizzare solo una coordinata invece di tutto il punto con le k coordinate??? RISPOSTA SI.
    struct kdtree_node *node = (struct kdtree_node *)malloc(sizeof(struct kdtree_node));
    node->median = ds[indexSorted[indexMedian] * k + cut];
    node->h_min = ds[indexSorted[start] * k + cut];   //valore di coordinata più piccola
    node->h_max = ds[indexSorted[end - 1] * k + cut]; //valore di coordinata più piccola

    // printf("\nPunto del nodo:  ");
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%f, ", ds[indexSorted[indexMedian] * k + i]);
    // }
    int numEleSx = indexMedian - start;
    int numEleDx = end - indexMedian - 1;
    // printf("\nsizeSX= %d  sizeDX= %d", numEleSx, numEleDx);

    if (numEleSx == 0 && numEleDx == 0)
    {
        // printf("\nNODO FOGLIA");
        node->left = NULL;
        node->right = NULL;
        return NULL;
    }
    else if (numEleSx == 0)
    {
        buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1);
        return node;
    }
    else if (numEleDx == 0)
    {
        node->right = NULL;
        buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0);
        return node;
    }
    else
    {
        buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0);
        buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1);
        return node;
    }
}

/*
*   buildTreeRoot server per costruire solo la radice del hdtree
*   end è il numero di elementi quindi va fatto end-1 per rpendere l'ultimo indice
*/
struct kdtree_node *buildTreeRoot(MATRIX ds, int *indexSorted, int liv, int end, int k)
{
    if (end <= 0)
    {
        printf("\nDATASET NULLO\n");
        return NULL;
    }
    int cut = liv % k; //variabile di cut per indice colonna da usare
    printf("\nCOSTRIUSCO RADICE\t livello= %d, size= %d, k= %d, cut= %d", liv, end, k, cut);

    quicksort(ds, indexSorted, cut, k, 0, end);

    // MyprintArray(ds, indexSorted, 0, size, k, cut);
    int indexMedian = findMedian(ds, indexSorted, 0, (end - 1) / 2, k, cut);

    printf("\nvalore MEDIANO= %f\t indexMedian=%d\t val puntoMin= %f\t val puntMax= %f ", ds[indexSorted[indexMedian] * k + cut], indexMedian, ds[indexSorted[0] * k + cut], ds[indexSorted[end - 1] * k + cut]);
    struct kdtree_node *root = (struct kdtree_node *)malloc(sizeof(struct kdtree_node));
    //SUGGERIMENTO: potrebbe bastare memorizzare solo una coordinata invece di tutto il punto con le k coordinate???
    // printf("\nPunto del nodo:  ");
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%f, ", ds[indexSorted[indexMedian] * k + i]);
    // }
    root->median = ds[indexSorted[indexMedian] * k + cut];
    root->h_min = ds[indexSorted[0] * k + cut];       //valore di coordinata più piccola
    root->h_max = ds[indexSorted[end - 1] * k + cut]; //valore di coordinata più grande

    int numEleSx = indexMedian;
    int numEleDx = end - indexMedian - 1;

    printf("\nsizeSX= %d  sizeDX= %d", numEleSx, numEleDx);
    root->left = buildTree(ds, indexSorted, liv + 1, 0, indexMedian, numEleSx, k, 0);
    root->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1);

    return root;
}

/*
*	K-d-Tree
* 	======================
*/
void kdtree(params *input)
{
    printf("\nInizio kdtree");
    // printf("\ndataset size%d, dataset k%d\n", input->n, input->k);
    int *indexSorted = (int *)malloc(input->n * sizeof(int)); //vettore che conterra indice riga dei punti ordinati
    struct kdtree_node *arrayTree = malloc(input->n * sizeof(struct kdtree_node));
    if (indexSorted == NULL || arrayTree == NULL)
    {
        printf("\nNO MEMORIA\n");
        exit(1);
    }
    for (int i = 0; i < input->n; i++)
    {
        indexSorted[i] = i;
    }
    input->kdtree = buildTreeRoot(input->ds, indexSorted, 0, input->n, input->k);
    //bisogna liberare la memoria
    free(indexSorted);
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

    printf("\n\ntime= %f seconds\n", time);
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
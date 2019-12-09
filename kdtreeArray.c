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
void MyprintArray(MATRIX, int *, int, int, int, int);
void quicksort(MATRIX, int *, int, int, int, int);
int findMedian(MATRIX, int *, int, int, int, int);
struct kdtree_node *buildTreeRoot(MATRIX, struct kdtree_node *, int, int *, int, int, int);
struct kdtree_node *buildTree(MATRIX, struct kdtree_node *, int, int *, int, int, int, int, int, int);

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

struct
{
    float median; //median coordinata
    float k_min, k_max;
    int indexVector, indexFatherVector, indexMedianPoint;
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
        // if (i == 100)
        //     break;
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
    for (int i = indexMedian - 1; i > start; i--)
    {
        if (dataset[indexSorted[i] * k + cut] < median)
        {
            return i + 1;
        }
    }
    return start;
}
int cont = 1;
/*
*   buildTree server per costruire tutti i nodi del kdtree
*   il valore end deve essere escluso indici vanno da [start, end)
*/
struct kdtree_node *buildTree(MATRIX ds, struct kdtree_node *arrayTree, int index, int *indexSorted, int liv, int start, int end, int numEle, int k, int type)
{

    int cut = liv % k; //variabile di cut per indice colonna da usare
    cont++;
    // type=0 nodo sinistro type=1 node destro
    if (type == 0)
        printf("\n\nCOSTRIUSCO FIGLIO SINISTRO liv= %d, start= %d end= %d, cut= %d, numEle= %d, cont= %d", liv, start, end, cut, numEle, cont);
    else
        printf("\n\nCOSTRIUSCO FIGLIO DESTRO liv= %d, start= %d end= %d, cut= %d, numEle= %d, cont= %d", liv, start, end, cut, numEle, cont);

    // if (cut == 5)
    // {
    //     // printf("\nFINE\n");
    //     return NULL;
    // }
    quicksort(ds, indexSorted, cut, k, start, end);

    int indexMedian = findMedian(ds, indexSorted, start, start + ((end - 1 - start) / 2), k, cut);

    // printf("\nvalore MEDIANO= %f\t indexMedian=%d\t val puntoMin= %f\t val puntMax= %f ", ds[indexSorted[indexMedian] * k + cut], indexMedian, ds[indexSorted[start] * k + cut], ds[indexSorted[end - 1] * k + cut]);

    // struct kdtree_node *node = (struct kdtree_node *)malloc(sizeof(struct kdtree_node));
    arrayTree[index]->median = ds[indexSorted[indexMedian] * k + cut];
    arrayTree[index]->k_min = ds[indexSorted[start] * k + cut];   //valore di coordinata più piccola
    arrayTree[index]->k_max = ds[indexSorted[end - 1] * k + cut]; //valore di coordinata più piccola
    arrayTree[index]->indexFatherVector = (index - 1) / 2;
    arrayTree[index]->indexMedianPoint = indexSorted[indexMedian];

    // printf("\nPunto del nodo:  ");
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%f, ", ds[indexSorted[indexMedian] * k + i]);
    // }
    int newSizeSx = indexMedian - start;
    int newSizeDx = end - indexMedian;
    // printf("\nsizeSX= %d  sizeDX= %d", newSizeSx, newSizeDx - 1);
    if (newSizeSx == 0 && newSizeDx - 1 == 0)
    {
        return NULL;
    }
    else if (newSizeSx == 0)
    {
        buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, indexMedian + 1, end, newSizeDx - 1, k, 1);
        return arrayTree;
    }
    else if (newSizeDx - 1 == 0)
    {
        buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, start, indexMedian, newSizeSx - 1, k, 0);
        return arrayTree;
    }
    else
    {
        buildTree(ds, arrayTree, (2 * index) + 1, indexSorted, liv + 1, start, indexMedian, newSizeSx, k, 0);
        buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, indexMedian + 1, end, newSizeDx - 1, k, 1);
        return arrayTree;
    }
}

/*
*   buildTreeRoot server per costruire solo la radice del hdtree
*   la size è il numero di elementi quindi va fatto size-1
*/
// struct kdtree_node *buildTreeRoot(MATRIX ds, int *indexSorted, int liv, int size, int k)
struct kdtree_node *buildTreeRoot(MATRIX ds, struct kdtree_node *arrayTree, int index, int *indexSorted, int liv, int end, int k)
{
    if (end <= 0)
    {
        printf("\nDATASET NULLO\n");
        return NULL;
    }
    int cut = liv % k; //variabile di cut per indice colonna da usare
    // printf("\nCOSTRIUSCO RADICE\t livello= %d, size= %d, k= %d, cut= %d", liv, size, k, cut);

    quicksort(ds, indexSorted, cut, k, 0, end);

    // MyprintArray(ds, indexSorted, 0, size, k, cut);
    int indexMedian = findMedian(ds, indexSorted, 0, (end - 1) / 2, k, cut);

    // printf("\nvalore MEDIANO= %f\t indexMedian=%d\t val puntoMin= %f\t val puntMax= %f ", ds[indexSorted[indexMedian] * k + cut], indexMedian, ds[indexSorted[0] * k + cut], ds[indexSorted[end - 1] * k + cut]);
    // struct kdtree_node *root = (struct kdtree_node *)malloc(sizeof(kdtree_node));
    // root->point = &ds[indexMedian]; //in point ci sarà l'indirizzo della prima cella del mediano
    // printf("\nPunto del nodo:  ");
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%f, ", ds[indexSorted[indexMedian] * k + i]);
    // }
    arrayTree[index]->median = ds[indexSorted[indexMedian] * k + cut];
    arrayTree[index]->k_min = ds[indexSorted[0] * k + cut];       //valore di coordinata più piccola
    arrayTree[index]->k_max = ds[indexSorted[end - 1] * k + cut]; //valore di coordinata più grande
    arrayTree[index]->indexFatherVector = index;
    arrayTree[index]->indexMedianPoint = indexSorted[indexMedian];

    int newSizeSx = indexMedian;
    int newSizeDx = end - newSizeSx;

    // printf("\nsizeSX= %d  sizeDX= %d", newSizeSx, newSizeDx);
    buildTree(ds, arrayTree, (2 * index) + 1, indexSorted, liv + 1, 0, indexMedian, newSizeSx, k, 0);
    buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, indexMedian + 1, end, newSizeDx - 1, k, 1);

    return arrayTree;
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
    struct kdtree_node *arrayTree = (struct kdtree_node *)malloc(input->n * sizeof(kdtree_node));
    if (indexSorted == NULL || arrayTree == NULL)
    {
        printf("\nNO MEMORIA\n");
        exit(1);
    }
    for (int i = 0; i < input->n; i++)
    {
        indexSorted[i] = i;
    }
    buildTreeRoot(input->ds, arrayTree, 0, indexSorted, 0, input->n, input->k);
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
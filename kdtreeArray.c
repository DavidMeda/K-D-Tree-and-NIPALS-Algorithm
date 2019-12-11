#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <malloc.h>

struct kdtree_node
{
    float median; //median coordinata
    float h_min, h_max;
    float *region;
    int indexMedianPoint;
};

#define MATRIX float *
#define KDTREE struct kdtree_node * // modificare con il tipo di dato utilizzato

void swap(int *, int *);
int partition(MATRIX, int *, int, int, int, int);
void MyprintArray(MATRIX, int *, int, int, int, int);
void quicksort(MATRIX, int *, int, int, int, int);
int findMedian(MATRIX, int *, int, int, int, int);
struct kdtree_node *buildTreeRoot(MATRIX, struct kdtree_node *, int, int *, int, int, int, float *, int);
struct kdtree_node *buildTree(MATRIX, struct kdtree_node *, int, int *, int, int, int, int, int, int, float *);
float *findRegion(MATRIX, int, int);

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
/*
    findMedian trova, per una stessa coordinata cut, il primo punto che abbia la stessa coordinata del mediano
    partendo dall'indice del mediano e risalendo fino al primo punto
*/
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

int cont = 1;
/*
*   buildTree server per costruire tutti i nodi del kdtree
*   il valore end deve essere escluso indici vanno da [start, end)
*   @type serve solo per stampare correttamente il figlio destro o sinistro
*/
struct kdtree_node *buildTree(MATRIX ds, struct kdtree_node *arrayTree, int index, int *indexSorted, int liv, int start, int end, int numEle, int k, int type, float *region)
{

    int cut = liv % k; //variabile di cut per indice colonna da usare
    cont++;
    // type=0 nodo sinistro type=1 node destro
    // if (type == 0)
    //     printf("\n\nCOSTRIUSCO FIGLIO SINISTRO liv= %d, start= %d end= %d, cut= %d, numEle= %d, cont= %d", liv, start, end, cut, numEle, cont);
    // else
    //     printf("\n\nCOSTRIUSCO FIGLIO DESTRO liv= %d, start= %d end= %d, cut= %d, numEle= %d, cont= %d", liv, start, end, cut, numEle, cont);

    arrayTree[index].h_min = ds[indexSorted[start] * k + cut - 1];      //valore di coordinata più piccola
    arrayTree[index].h_max = ds[indexSorted[end - 1] * k + cut - 1];    //valore di coordinata più piccola
    region[2 * (cut - 1)] = ds[indexSorted[start] * k + cut - 1];       //aggiorno il valore di coordinata nell'intera regione
    region[2 * (cut - 1) + 1] = ds[indexSorted[end - 1] * k + cut - 1]; //aggiorno il valore di coordinata nell'intera regione
    arrayTree[index].region = region;
    
    // printf("  puntoMin= %f  val puntMax= %f ", ds[indexSorted[start] * k + cut - 1], ds[indexSorted[end - 1] * k + cut - 1]);
    
    //server per in debug
    // if (cut == 2)
    // {
    //     // printf("\nFINE\n");
    //     return NULL;
    // }

    // for (int z = 0; z < 6; z++)
    // {
    //     printf("\nmin= %f - max= %f  ", arrayTree[index].region[2 * z], arrayTree[index].region[(2 * z) + 1]);
    // }
    quicksort(ds, indexSorted, cut, k, start, end);
    int indexMedian = findMedian(ds, indexSorted, start, start + ((end - 1 - start) / 2), k, cut);
    // printf("\nvalore MEDIANO= %f\t indexMedian=%d\t", ds[indexSorted[indexMedian] * k + cut], indexMedian);

    arrayTree[index].median = ds[indexSorted[indexMedian] * k + cut]; //aggiungo il valore del mediano per la coordinata k (cut) nel nodo
    arrayTree[index].indexMedianPoint = indexSorted[indexMedian];     //memorizzo l'indice del dataset corrispondente al punto mediano inserito nel nodo
    // TODO: conviene usare l'indirizzo di memoria oppure l'indice della posizione nel dataset ??

    //stampa le coordinate del punto mediano
    // printf("\nPunto del nodo:  ");
    // for (int i = 0; i < k; i++)
    // {
    //     printf("%f, ", ds[indexSorted[indexMedian] * k + i]);
    // }
    int numEleSx = indexMedian - start;
    int numEleDx = end - indexMedian - 1;
    // printf("\nsizeSX= %d  sizeDX= %d", newSizeSx, newSizeDx);
    if (numEleSx == 0 && numEleDx == 0)
    {
        return NULL;
    }
    else if (numEleSx == 0)
    {
        buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1, region);
        return arrayTree;
    }
    else if (numEleDx == 0)
    {
        buildTree(ds, arrayTree, (2 * index) + 1, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0, region);
        return arrayTree;
    }
    else
    {
        buildTree(ds, arrayTree, (2 * index) + 1, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0, region);
        buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1, region);
        return arrayTree;
    }
}

/*
*   buildTreeRoot server per costruire solo la radice del hdtree
*   la size è il numero di elementi quindi va fatto size-1
*/
// struct kdtree_node *buildTreeRoot(MATRIX ds, int *indexSorted, int liv, int size, int k)
struct kdtree_node *buildTreeRoot(MATRIX ds, struct kdtree_node *arrayTree, int index, int *indexSorted, int liv, int end, int k, float *region, int indexRegion)
{
    if (end <= 0)
    {
        printf("\nDATASET NULLO\n");
        return NULL;
    }
    int cut = liv % k; //variabile di cut per indice colonna da usare
    // printf("\nCOSTRIUSCO RADICE\t livello= %d, size= %d, k= %d, cut= %d", liv, end, k, cut);

    quicksort(ds, indexSorted, cut, k, 0, end);

    // MyprintArray(ds, indexSorted, 0, size, k, cut);
    int indexMedian = findMedian(ds, indexSorted, 0, (end - 1) / 2, k, cut);

    // printf("\nvalore MEDIANO= %f\t indexMedian=%d\t val puntoMin= %f\t val puntMax= %f ", ds[indexSorted[indexMedian] * k + cut], indexMedian, ds[indexSorted[0] * k + cut], ds[indexSorted[end - 1] * k + cut]);
    // printf("\nPunto del nodo:  ");
    // for (int i = 0; i < k; i++)
    // {
    //     printf("%f, ", ds[indexSorted[indexMedian] * k + i]);
    // }
    arrayTree[index].median = ds[indexSorted[indexMedian] * k + cut];
    arrayTree[index].h_min = ds[indexSorted[0] * k + cut];       //valore di coordinata più piccola
    arrayTree[index].h_max = ds[indexSorted[end - 1] * k + cut]; //valore di coordinata più grande
    arrayTree[index].indexMedianPoint = indexSorted[indexMedian];

    region[2 * cut] = ds[indexSorted[0] * k + cut];
    region[(2 * cut) + 1] = ds[indexSorted[end - 1] * k + cut];
    arrayTree[index].region = region;
    // for (int z = 0; z < 7; z++)
    // {
    //     printf("\nmin= %f - max= %f  ", region[2 * z], region[(2 * z) + 1]);
    // }

    int numEleSx = indexMedian;
    int numEleDx = end - indexMedian - 1;

    // printf("\nsizeSX= %d  sizeDX= %d", numEleSx, numEleDx);
    buildTree(ds, arrayTree, (2 * index) + 1, indexSorted, liv + 1, 0, indexMedian, numEleSx, k, 0, region);
    buildTree(ds, arrayTree, (2 * index) + 2, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1, region);

    return arrayTree;
}

float *findRegion(MATRIX ds, int n, int k)
{
    float *region = (float *)malloc(k * 2 * sizeof(float));

    for (int j = 0; j < k; j++)
    {
        float h_min = ds[j], h_max = ds[j];
        for (int i = 0; i < n; i++)
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

/*
*	K-d-Tree
* 	======================
*/
void kdtree(params *input)
{
    printf("\nInizio kdtree");
    int *indexSorted = (int *)malloc(input->n * sizeof(int)); //vettore che conterra indice riga dei punti ordinati
    struct kdtree_node *arrayTree = (struct kdtree_node *)malloc(input->n * input->n * sizeof(struct kdtree_node));

    if (indexSorted == NULL || arrayTree == NULL)
    {
        printf("\nNO MEMORIA\n");
        exit(1);
    }
    for (int i = 0; i < input->n; i++)
    {
        indexSorted[i] = i;
    }
    float *region = findRegion(input->ds, input->n, input->k);

    input->kdtree = buildTreeRoot(input->ds, arrayTree, 0, indexSorted, 0, input->n, input->k, region, 0);

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
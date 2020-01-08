#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <malloc.h>
#include <libgen.h>

#define MATRIX float *
#define KDTREE struct kdtree_node * // modificare con il tipo di dato utilizzato

void swap(int *, int *);
int partition(MATRIX, int *, int, int, int, int);
void MyprintArray(MATRIX, int *, int, int, int, int);
void quicksort(MATRIX, int *, int, int, int, int);
int findMedian(MATRIX, int *, int, int, int, int);
struct kdtree_node *buildTreeRoot(MATRIX, int *, int, int, int);
struct kdtree_node *buildTree(MATRIX, int *, int, int, int, int, int, int);
float *findRegion(MATRIX, int, int);
float euclideanDistance(float *, MATRIX, int, int);
float euclideanDistanceDataset(MATRIX, int, MATRIX, int, int);
float distanceRoot(float *, MATRIX, int, int, float *);
float distanceChild(KDTREE, MATRIX, int, int, float *, int);
int *rangeQueryChild(MATRIX, KDTREE, int, int, MATRIX, int, int, int, float *, int *);
int *rangeQueryRoot(MATRIX, KDTREE, int, int, MATRIX, int, int, int, float *, float *, int *);
void centraMatrice(MATRIX, int, int);
float calcolaT(float *, int, int, int);
void prodottoMatrice(float *, int, int, MATRIX, int, int, float *, int, int, int, int, int);
void prodottoMatriceTrasp(float *, int, int, MATRIX, int, int, float *, int, int, int, int, int);
float norma(float *, int, int, int);
void dividi(float *, int, int, int, float);
void aggiornaDataset(MATRIX, int, int, float *, int, int, float *, int, int, int, int, int);

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
    float h_min, h_max;
    int indexMedianPoint;
    int numPoint; //forse non serve

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

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
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
        fwrite(&k, 4, 1, fp);
        fwrite(&n, 4, 1, fp);
        for (i = 0; i < n; i++)
        {
            fwrite(X, 4, k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += 4 * k;
        }
    }
    fclose(fp);
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

float euclideanDistance(float *point, MATRIX q, int indexQ, int k)
{

    float sum = 0;
    int i;
    for (i = 0; i < k; i++)
    {
        sum += (point[i] - q[indexQ * k + i]) * (point[i] - q[indexQ * k + i]);
    }
    return sqrtf(sum);
}

float euclideanDistanceDataset(MATRIX ds, int indexMedian, MATRIX q, int indexQ, int k)
{

    float sum = 0;
    int i;
    for (i = 0; i < k; i++)
    {
        sum += (ds[indexMedian * k + i] - q[indexQ * k + i]) * (ds[indexMedian * k + i] - q[indexQ * k + i]);
    }
    return sqrtf(sum);
}

//distanza tra il punto q del querySet e l'intera regione indicizzata
float distanceRoot(float *region, MATRIX q, int indexQ, int k, float *point)
{
    int j;
    for (j = 0; j < k; j++)
    {
        if (q[indexQ * k + j] <= region[2 * j])
            point[j] = region[2 * j];
        else if (q[indexQ * k + j] >= region[(2 * j) + 1])
            point[j] = region[(2 * j) + 1];
        else
            point[j] = q[indexQ * k + j];
    }

    return euclideanDistance(point, q, indexQ, k);
}

float distanceChild(KDTREE node, MATRIX q, int indexQ, int k, float *point, int cut)
{
    if (q[indexQ * k + cut] <= node->h_min)
        point[cut] = q[indexQ * k + cut];
    else if (q[indexQ * k + cut] >= node->h_max)
        point[cut] = q[indexQ * k + cut];
    else
        point[cut] = q[indexQ * k + cut];

    return euclideanDistance(point, q, indexQ, k);
}

/*
*   buildTree serve per costruire tutti i nodi del kdtree
*   il valore end deve essere escluso indici vanno da [start, end) quindi end è escluso l'ultimo elemento si trova a end-1
*/
struct kdtree_node *buildTree(MATRIX ds, int *indexSorted, int liv, int start, int end, int numEle, int k, int type)
{
    if (numEle == 0)
        return NULL;

    int cut = liv % k; //variabile di cut per indice colonna da usare
    struct kdtree_node *node = (struct kdtree_node *)get_block(sizeof(struct kdtree_node), 1);
    node->h_min = ds[indexSorted[start] * k + (cut - 1)];   //valore di coordinata più piccola per il padre
    node->h_max = ds[indexSorted[end - 1] * k + (cut - 1)]; //valore di coordinata più piccola per il padre
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
        node->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1);
        return node;
    }
    else if (numEleDx == 0)
    {
        node->right = NULL;
        node->left = buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0);
        return node;
    }
    else
    {
        node->left = buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, 0);
        node->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1);
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

    quicksort(ds, indexSorted, cut, k, 0, end);

    int indexMedian = findMedian(ds, indexSorted, 0, (end - 1) / 2, k, cut);

    struct kdtree_node *root = (struct kdtree_node *)get_block(sizeof(struct kdtree_node), 1);

    root->medianCoordinate = ds[indexSorted[indexMedian] * k + cut]; //valore di coordinata del punto mediano
    root->indexMedianPoint = indexSorted[indexMedian];               //indice del punto mediano nel dataset
    root->h_min = ds[indexSorted[0] * k + cut];                      //valore di coordinata più piccola
    root->h_max = ds[indexSorted[end - 1] * k + cut];                //valore di coordinata più grande
    root->numPoint = end;

    int numEleSx = indexMedian;
    int numEleDx = end - indexMedian - 1;

    root->left = buildTree(ds, indexSorted, liv + 1, 0, indexMedian, numEleSx, k, 0);
    root->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, 1);

    return root;
}

float *findRegion(MATRIX ds, int n, int k)
{
    float *region = (float *)get_block(sizeof(float), k * 2);
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

/*  Questo metodo effettua la rangeQuery sui sottoalberi figli di root, la distanza della region viene effettuata solo per le coordinate di cut del nodo padre
*   poichè le restanti coordinate non cambiano e alla fine si eseguono le chiamate ricorsive sui figli finchè non si arriva a un nodo foglia
*/
int *rangeQueryChild(MATRIX ds, KDTREE node, int n, int k, MATRIX q, int indexQ, int r, int liv, float *point, int *list)
{
    int cut = (liv - 1) % k; //il cut del nodo padre

    if (distanceChild(node, q, indexQ, k, point, cut) > r)
    {
        printf("stop  ");
        return list;
    }

    if (euclideanDistanceDataset(ds, node->indexMedianPoint, q, indexQ, k) <= r)
    {
        list[indexQ * n + indexList] = node->indexMedianPoint;

        indexList++;
    }

    if (node->left != NULL)
    {
        rangeQueryChild(ds, node->left, n, k, q, indexQ, r, liv + 1, point, list);
    }
    if (node->right != NULL)
    {
        rangeQueryChild(ds, node->right, n, k, q, indexQ, r, liv + 1, point, list);
    }
    return list;
}

/*  Questo metodo server per effettuare il controllo tra il punto di query e la regione indicizzata dall'intero dataset
*   viene usato region che contiene il vettore delle k coppia di (h_min, h_max)
*   alla fine il metodo richiama la ricerca sui suoi figli
*/
int *rangeQueryRoot(MATRIX ds, KDTREE root, int n, int k, MATRIX q, int indexQ, int r, int liv, float *region, float *point, int *list)
{

    if (distanceRoot(region, q, indexQ, k, point) > r)
    {
        return NULL;
    }

    if (euclideanDistanceDataset(ds, root->indexMedianPoint, q, indexQ, k) <= r)
    {
        list[indexQ * n + indexList] = root->indexMedianPoint;
        indexList++;
    }

    if (root->left != NULL)
    {
        rangeQueryChild(ds, root->left, n, k, q, indexQ, r, liv + 1, point, list);
    }
    if (root->right != NULL)
    {
        rangeQueryChild(ds, root->right, n, k, q, indexQ, r, liv + 1, point, list);
    }
    return list;
}

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


void prodottoMatriceTrasp(float *result, int rigaRes, int cut, MATRIX ds, int rigaA, int colA, float *vect, int rigaB, int colB, int k, int n, int h)
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
                sum += ds[j * k + m] * vect[j * h + cut];
            }
            result[m * h + cut] = sum;
        }
    }
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
            }
            result[m * h + cut] = sum;
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

    int m, i, j;
    float sum = 0;
    for (m = 0; m < rigaA; m++)
    {
        sum = 0;
        for (i = 0; i < colB; i++)
        {
            sum += u[m * h + cut] * v[i * h + cut];
        }
        ds[m * k + i] -= sum;
    }

}

/*
*	PCA
* 	=====================
*/
void pca(params *input)
{
    printf("\nINIZIO PCA");

    centraMatrice(input->ds, input->n, input->k);
    float theta = 1 * exp(-8);
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
    for (i = 0; i < input->h; i++)
    {
        do
        {
            prodottoMatriceTrasp(input->V, input->k, i, input->ds, input->k, input->n, input->U, input->n, 1, input->k, input->n, input->h);

            t = calcolaT(input->U, input->n, input->h, i);
            dividi(input->V, input->k, input->h, i, t);
            norm = norma(input->V, input->k, input->h, i);
            dividi(input->V, input->k, input->h, i, norm);

            prodottoMatrice(input->U, input->n, i, input->ds, input->n, input->k, input->V, input->k, 1, input->k, input->n, input->h);
            tempV = calcolaT(input->V, input->k, input->h, i);
            dividi(input->U, input->n, input->h, i, tempV);
            t1 = calcolaT(input->U, input->n, input->h, i);

        } while (t1 - t >= theta * t1);

        aggiornaDataset(input->ds, input->n, input->k, input->U, input->n, 1, input->V, 1, input->k, input->h, input->k, i);
    }
    printf("\nFINE PCA");
}

/*
*	K-d-Tree
* 	======================
*/
void kdtree(params *input)
{

    printf("\nInizio kdtree");
    // printf("\ndataset size%d, dataset k%d\n", input->n, input->k);
    int *indexSorted = (int *)get_block(sizeof(int), input->n); //vettore che conterra indice riga dei punti ordinati

    if (indexSorted == NULL)
    {
        printf("\nNO MEMORIA\n");
        exit(-2);
    }
    int i;
    for (i = 0; i < input->n; i++)
    {
        indexSorted[i] = i;
    }
    // si può usare memset
    // memset(indexSorted,0,input->n*sizeof(int));

    input->region = findRegion(input->ds, input->n, input->k);
    input->kdtree = buildTreeRoot(input->ds, indexSorted, 0, input->n, input->k);

    // // printTree(input->kdtree);

    //bisogna liberare la memoria

    free_block(indexSorted);
    printf("\nfine KDTREE");

    // free(region);
}

void range_query(params *input)
{
    printf("\n inizio QUERY");

    float *point = (float *)get_block(sizeof(float), input->k); //punto che conterrà il punto più vicino alla region per il confronto
    input->QA = (int *)get_block(sizeof(int), input->n * input->nq);
    if (point == NULL || input->QA == NULL)
    {
        printf("\nNO MEMORIA");
        exit(-2);
    }
    int i;
    for (i = 0; i < input->nq; i++)
    {
        rangeQueryRoot(input->ds, input->kdtree, input->n, input->k, input->qs, i, input->r, 0, input->region, point, input->QA);
        input->QA[i * input->n + indexList] = -1;
        input->nQA += indexList + 1;
        indexList = 0;
    }
    free_block(point);
}

int main(int argc, char const *argv[])
{

    char fname[256];
    char *dsname;
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
    input->display = 0;
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

    dsname = basename(strdup(input->filename));

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

    if(input->h > 0){
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
        printf("\n %s \n",dsname);
        sprintf(fname, "%s.U", dsname);
        sprintf(fname, "%s.V", dsname);
    }else
        time = -1;
       
    if (!input->silent)
        printf("\nPCA time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    
    //
    // Costruzione K-d-Tree
    //
    
    if(input->kdtree_enabled){
        t = clock();
        kdtree(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
    }else
        time = -1;
    if (!input->silent)
        printf("\nIndexing time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Range query search
    //
    
    if(input->r >= 0){
        t = clock();
        range_query(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
    }else
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
            int flag = 0;
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
        save_data(fname, input->QA, input->nQA, 2);
    }

    if (!input->silent)
        printf("\nDone.\n");

    return 0;
}
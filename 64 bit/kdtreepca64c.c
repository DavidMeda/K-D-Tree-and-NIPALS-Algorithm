#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define MATRIX float *
#define KDTREE struct kdtree_node * // modificare con il tipo di dato utilizzato

void swap(int *, int *);
int partition(MATRIX, int *, int, int, int, int);
void quicksort(MATRIX, int *, int, int, int, int);
int findMedian(MATRIX, int *, int, int, int, int);
struct kdtree_node *buildTreeRoot(MATRIX, int *, int, int, int);
struct kdtree_node *buildTree(MATRIX, int *, int, int, int, int, int, float *);
float *findRegion(MATRIX, int, int);
float euclideanDistance(MATRIX qs, int id_qs, MATRIX ds, int id_ds, int k);
int rangeQuery(MATRIX, struct kdtree_node *, MATRIX, int, float, int, int, int *, float *, int *);
void centraMatrice(MATRIX, int, int, MATRIX, int);
float calcolaT(float *vect, int numEle);
void prodottoMatrice(float *u, MATRIX ds, int rigaDS, float *v, int k);
void prodottoMatriceTrasp(float *v, MATRIX ds, float *u, int numEleU, int k);
float norma(float *vect, int numEle);
void dividi(float *vect, int numEle, float value);
void aggiornaDataset(MATRIX ds, int n, int k, float *u, float *v);
float *calcoloQ(MATRIX, MATRIX, int, int, int);
int indexList = 0, kOld = 0;

// funzioni Assembly  AVX
extern float euclideanDistanceAss_64(MATRIX ds, MATRIX qs, int k);
extern void dividiAss_64(float *vett, int numEle, float *value);
extern void aggiornaDatasetAss_64(MATRIX ds, float *u, float *v, int n, int k);
extern void prodMatriceAss_64(float *ds, float *u, float *v, int n, int k);
extern void prodMatriceTrasAss_64(float *ds, float *u, float *v, int n, int k);

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
    int **QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA;  //numero di risposte alle query
} params;

struct kdtree_node
{
    float medianCoordinate; //valore di coordinata cut per il punto mediano
    int indexMedianPoint;   //indice di riga del punto mediano
    float *region;          //regione indicizzata da root formata da k coppie (h_min, h_max)
    struct kdtree_node *left, *right;
};

void *get_block(int size, int elements)
{
    return _mm_malloc(elements * size, 32);
}

void free_block(void *p)
{
    _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols)
{
    return (MATRIX)get_block(sizeof(float), rows * cols);
}

void dealloc_matrix(MATRIX mat)
{
    free_block(mat);
}

void deallocate(KDTREE tree)
{
    if (tree == NULL)
        return;
    deallocate(tree->right);
    deallocate(tree->left);
    free_block(tree->region);
    free_block(tree);
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

void saveDataVicini(char *filename, int **X, int puntiTot, int n_qs, int n_ds)
{
    FILE *fp;
    int i;
    fp = fopen(filename, "wb");
    int temp[2];
    if (X != NULL)
    {
        int col = 2;
        fwrite(&col, sizeof(int), 1, fp);
        fwrite(&puntiTot, sizeof(int), 1, fp);

        for (i = 0; i < n_qs; i++)
        {
            if (X[i] != NULL)
            {
                temp[0] = i;
                for (int j = 0; j < n_ds; j++)
                {
                    if (X[i][j] == -1)
                        break;
                    temp[1] = X[i][j];
                    fwrite(temp, sizeof(temp), 1, fp);
                }
            }
        }
    }
    fclose(fp);
}

/*  metodi KDTREE
*
*/

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

struct kdtree_node *buildTree(MATRIX ds, int *indexSorted, int liv, int start, int end, int numEle, int k, float *regionp)
{
    if (numEle == 0)
        return NULL;

    int cut = liv % k; //variabile di cut per indice colonna da usare
    struct kdtree_node *node = (struct kdtree_node *)get_block(sizeof(struct kdtree_node), 1);
    node->region = get_block(sizeof(float), 2 * k);
    if (node == NULL || node->region == NULL)
    {
        printf("NO MEMORIA NODO\n");
        exit(1);
    }
    memcpy(node->region, regionp, sizeof(float) * 2 * k);
    int oldCut = (liv - 1) % k;
    node->region[2 * (oldCut)] = ds[indexSorted[start] * k + (oldCut)];
    node->region[2 * (oldCut) + 1] = ds[indexSorted[end - 1] * k + (oldCut)];

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
        node->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, node->region);
        return node;
    }
    else if (numEleDx == 0)
    {
        node->right = NULL;
        node->left = buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, node->region);
        return node;
    }
    else
    {
        node->left = buildTree(ds, indexSorted, liv + 1, start, indexMedian, numEleSx, k, node->region);
        node->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, node->region);
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

    struct kdtree_node *root = (struct kdtree_node *)get_block(sizeof(struct kdtree_node), 1);
    if (root == NULL)
    {
        printf("NO MEMORIA NODO\n");
        exit(1);
    }
    root->medianCoordinate = ds[indexSorted[indexMedian] * k + cut]; //valore di coordinata del punto mediano
    root->indexMedianPoint = indexSorted[indexMedian];               //indice del punto mediano nel dataset
    root->region = findRegion(ds, end, k);
    int numEleSx = indexMedian;
    int numEleDx = end - indexMedian - 1;

    root->left = buildTree(ds, indexSorted, liv + 1, 0, indexMedian, numEleSx, k, root->region);
    root->right = buildTree(ds, indexSorted, liv + 1, indexMedian + 1, end, numEleDx, k, root->region);

    return root;
}

float *findRegion(MATRIX ds, int n, int k)
{
    float *region = (float *)get_block(sizeof(float), 2 * k);
    float h_min = 0, h_max = 0;
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

////////FINE KDTREE METHOD

/*  METODI PCA
*
*/

void moltiplica(MATRIX q, MATRIX V, float *c, int nq, int k, int h)
{
    int i, j, z;
    for (i = 0; i < nq; i++)
    {
        for (j = 0; j < h; j++)
        {
            for (z = 0; z < k; z++)
            {
                c[i * h + j] += q[i * k + z] * V[z * h + j];
            }
        }
    }
}

void centraMatrice(MATRIX ds, int n, int k, MATRIX qs, int nq)
{
    int i, j;
    float acc = 0;
    float mean = 0;
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
        for (int z = 0; z < nq; z++)
        {
            qs[z * k + j] = qs[z * k + j] - mean;
        }
    }
}

float calcolaT(float *vect, int numEle)
{
    float res = 0;
    for (int i = 0; i < numEle; i++)
    {
        res += (vect[i]) * (vect[i]);
    }
    return res;
}

float norma(float *vect, int numEle)
{
    float acc = 0;
    for (int i = 0; i < numEle; i++)
    {
        acc += vect[i] * vect[i];
    }
    return sqrtf(acc);
}

void prodottoMatriceTrasp(float *v, MATRIX ds, float *u, int numEleU, int k)
{
    // memset(v, 0, sizeof(float) * k);
    // prodMatriceTrasAss_64(ds, u, v, numEleU, k);

    // int i, j;
    // float sum = 0;
    // for (i = 0; i < k; i++)
    // {
    //     sum = 0;
    //     for (j = 0; j < numEleU; j++)
    //     {
    //         sum += ds[j * k + i] * u[j];
    //     }
    //     v[i] = sum;
    // }
}

void prodottoMatrice(float *u, MATRIX ds, int rigaDS, float *v, int k)
{
    // prodMatriceAss_64(ds, u, v, rigaDS, k);

    // int i, j;
    // float sum = 0;
    // for (i = 0; i < rigaDS; i++)
    // {
    //     sum = 0;
    //     for (j = 0; j < k; j++)
    //     {
    //         sum += ds[i * k + j] * v[j];
    //     }
    //     u[i] = sum;
    // }
}

void dividi(float *vect, int numEle, float value)
{

    // for (int i = 0; i < numEle; i++)
    // {
    //     vect[i] = vect[i] / value;
    // }
}

void aggiornaDataset(MATRIX ds, int n, int k, float *u, float *v)
{
    // aggiornaDatasetAss_64(ds, u, v, n, k);

    // int i, j;
    // for (i = 0; i < n; i++)
    // {
    //     for (j = 0; j < k; j++)
    //     {
    //         ds[i * k + j] -= u[i] * v[j];
    //     }
    // }
}

float *calcoloQ(MATRIX q, MATRIX V, int nq, int k, int h)
{
    float *q1 = (float *)get_block(sizeof(float), h * nq);
    if (q1 == NULL)
    {
        printf("NO MEMORIA q1\n");
        exit(1);
    }
    moltiplica(q, V, q1, nq, k, h);
    return q1;
}

/////FINE METODI PCA

void pca(params *input)
{
    centraMatrice(input->ds, input->n, input->k, input->qs, input->nq);
    input->V = get_block(sizeof(float), input->k * input->h); //dimensioni (k x h)
    input->U = get_block(sizeof(float), input->n * input->h); // dimensioni (n x h)
    float *u = get_block(sizeof(float), input->n);
    float *v = get_block(sizeof(float), input->k);

    if (input->V == NULL || input->U == NULL || u == NULL || v == NULL)
    {
        printf("\nNo MEMORIA PCA\n");
        exit(1);
    }
    int i;
    for (i = 0; i < input->n; i++)
    {
        u[i] = input->ds[i * input->k];
    }
    float theta = 1e-8;
    float norm = 0, tempV = 0, diff = 0, t = 0, t1 = 0;
    for (i = 0; i < input->h; i++)
    {

        diff = 0, t = 0, t1 = 0;

        do
        {
            memset(v, 0, sizeof(float) * input->k);
            prodMatriceTrasAss_64(input->ds, u, v, input->n, input->k);
            // prodottoMatriceTrasp(v, input->ds, u, input->n, input->k);

            t = calcolaT(u, input->n);
            // dividi(v, input->k, t);
            dividiAss_64(v, input->k, &t);
            norm = norma(v, input->k);

            // dividi(v, input->k, norm);
            dividiAss_64(v, input->k, &norm);

            prodMatriceAss_64(input->ds, u, v, input->n, input->k);

            // prodottoMatrice(u, input->ds, input->n, v, input->k);
            tempV = calcolaT(v, input->k);
            // dividi(u, input->n, tempV);
            dividiAss_64(u, input->n, &tempV);

            t1 = calcolaT(u, input->n);

            diff = t1 - t;
            if (diff < 0)
                diff = diff * -1;

        } while (diff >= theta * t1);

        // aggiornaDataset(input->ds, input->n, input->k, u, v);
        aggiornaDatasetAss_64(input->ds, u, v, input->n, input->k);

        for (int j = 0; j < input->k; j++)
        {
            input->V[j * input->h + i] = v[j];
        }

        for (int z = 0; z < input->n; z++)
        {
            input->U[z * input->h + i] = u[z];
        }
    }

    free_block(u);
    free_block(v);
    free_block(input->ds);
    input->ds = input->U;
    float *newQS = calcoloQ(input->qs, input->V, input->nq, input->k, input->h);
    kOld = input->k;
    input->k = input->h;
    free_block(input->qs);
    input->qs = newQS;
}

void kdtree(params *input)
{

    int *indexSorted = get_block(sizeof(int), input->n);

    if (indexSorted == NULL)
    {
        printf("\nNO MEMORIA indexSorted\n");
        exit(1);
    }
    int i;
    for (i = 0; i < input->n; i++)
    {
        indexSorted[i] = i;
    }

    input->kdtree = buildTreeRoot(input->ds, indexSorted, 0, input->n, input->k);

    free_block(indexSorted);
}

/*  METODI RANGE QUERY
*
*/

float distance(float *h, MATRIX qs, int id_qs, int k, float *p)
{

    int repeat = k / 4, left = k % 4, i = 0;

    while (repeat--)
    {
        //blocco 1
        if (qs[k * id_qs + i] <= h[i * 2])
            p[i] = h[i * 2];
        else if (qs[k * id_qs + i] >= h[(i * 2) + 1])
            p[i] = h[(i * 2) + 1];
        else
            p[i] = qs[k * id_qs + i];

        //blocco2
        if (qs[k * id_qs + (i + 1)] <= h[(i + 1) * 2])
            p[(i + 1)] = h[(i + 1) * 2];
        else if (qs[k * id_qs + (i + 1)] >= h[((i + 1) * 2) + 1])
            p[(i + 1)] = h[((i + 1) * 2) + 1];
        else
            p[(i + 1)] = qs[k * id_qs + (i + 1)];

        //blocco 3
        if (qs[k * id_qs + (i + 2)] <= h[(i + 2) * 2])
            p[(i + 2)] = h[(i + 2) * 2];
        else if (qs[k * id_qs + (i + 2)] >= h[((i + 2) * 2) + 1])
            p[(i + 2)] = h[((i + 2) * 2) + 1];
        else
            p[(i + 2)] = qs[k * id_qs + (i + 2)];

        //blocco 4
        if (qs[k * id_qs + (i + 3)] <= h[(i + 3) * 2])
            p[(i + 3)] = h[(i + 3) * 2];
        else if (qs[k * id_qs + (i + 3)] >= h[((i + 3) * 2) + 1])
            p[(i + 3)] = h[((i + 3) * 2) + 1];
        else
            p[(i + 3)] = qs[k * id_qs + (i + 3)];

        i += 4;
    }

    switch (left)
    {
    case 3:
        if (qs[k * id_qs + (i + 2)] <= h[(i + 2) * 2])
            p[(i + 2)] = h[(i + 2) * 2];
        else if (qs[k * id_qs + (i + 2)] >= h[((i + 2) * 2) + 1])
            p[(i + 2)] = h[((i + 2) * 2) + 1];
        else
            p[(i + 2)] = qs[k * id_qs + (i + 2)];

    case 2:
        if (qs[k * id_qs + (i + 1)] <= h[(i + 1) * 2])
            p[(i + 1)] = h[(i + 1) * 2];
        else if (qs[k * id_qs + (i + 1)] >= h[((i + 1) * 2) + 1])
            p[(i + 1)] = h[((i + 1) * 2) + 1];
        else
            p[(i + 1)] = qs[k * id_qs + (i + 1)];

    case 1:
        if (qs[k * id_qs + i] <= h[i * 2])
            p[i] = h[i * 2];
        else if (qs[k * id_qs + i] >= h[(i * 2) + 1])
            p[i] = h[(i * 2) + 1];
        else
            p[i] = qs[k * id_qs + i];
    case 0:;
    }
    return euclideanDistanceAss_64(&qs[id_qs * k], p, k);
}

int rangeQuery(MATRIX ds, struct kdtree_node *tree, MATRIX qs, int id_qs, float r, int k, int n, int *list, float *point, int *nQA)
{

    if (tree == NULL || distance(tree->region, qs, id_qs, k, point) > r)
    {
        return 0;
    }

    if (euclideanDistanceAss_64(&qs[id_qs * k], &ds[tree->indexMedianPoint * k], k) <= r)
    {
        list[indexList] = tree->indexMedianPoint;
        indexList++;
        *nQA = *nQA + 1;
    }
    if (tree->left != NULL)
    {
        rangeQuery(ds, tree->left, qs, id_qs, r, k, n, list, point, nQA);
    }
    if (tree->right != NULL)
    {
        rangeQuery(ds, tree->right, qs, id_qs, r, k, n, list, point, nQA);
    }

    return 0;
}

////FINE METODI RANGE QUERY

void range_query(params *input)
{
    input->QA = get_block(sizeof(int *), input->nq);
    float *point = get_block(sizeof(float), input->k);
    if (input->QA == NULL)
    {
        printf("\nNO MEMORIA QA\n");
        exit(1);
    }
    if (point == NULL)
    {
        printf("\nNO MEMORIA point\n");
        exit(1);
    }
    int i;
    for (i = 0; i < input->nq; i++)
    {
        input->QA[i] = get_block(sizeof(int), input->n);
        rangeQuery(input->ds, input->kdtree, input->qs, i, input->r, input->k, input->n, input->QA[i], point, &input->nQA);
        if (indexList == 0)
        {
            free_block(input->QA[i]);
            input->QA[i] = NULL;
        }
        else
        {
            input->QA[i][indexList] = -1;
        }
        indexList = 0;
    }
    free_block(point);
}

int main(int argc, char **argv)
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

    if (input->h > 0)
    {
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t) / CLOCKS_PER_SEC;
        sprintf(fname, "%s.U", dsname);
        //salvo U di dimensioni (n, h)
        save_data(fname, input->U, input->n, input->k);
        sprintf(fname, "%s.V", dsname);
        //salvo V di dimensioni (k, h)
        save_data(fname, input->V, kOld, input->k);
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

    if (input->kdtree_enabled)
    {
        t = clock();
        kdtree(input);
        t = clock() - t;
        time = ((float)t) / CLOCKS_PER_SEC;
    }
    else
        time = -1;
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
            for (i = 0; i < input->nq; i++)
            {
                if (input->QA[i] != NULL)
                {
                    for (j = 0; j < input->n; j++)
                    {
                        if (input->QA[i][j] == -1)
                            break;
                        else
                        {
                            printf("%d,  ", i);
                            printf("%d\n", input->QA[i][j]);
                        }
                    }
                }
            }
        }
        // printf("\n %d ", input->nQA);
        sprintf(fname, "%s.qa", dsname);
        saveDataVicini(fname, input->QA, input->nQA, input->nq, input->n);
    }

    if (!input->silent)
        printf("\nDone.\n");

    if (input->QA != NULL)
    {
        for (int i = 0; i < input->nq; i++)
        {
            if (input->QA[i] != NULL)
            {
                free_block(input->QA[i]);
            }
        }
        free_block(input->QA);
    }
    if (input->U != NULL)
    {
        free_block(input->U);
    }
    if (input->V != NULL)
    {
        free_block(input->V);
    }
    if (input->ds != NULL && input->h <= 0)
    {
        free_block(input->ds);
    }
    if (input->qs != NULL)
    {
        free_block(input->qs);
    }

    if (input->kdtree != NULL)
    {
        deallocate(input->kdtree);
    }

    free(input);

    return 0;
}
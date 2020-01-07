#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <malloc.h>

#define MATRIX float *
#define KDTREE struct kdtree_node * // modificare con il tipo di dato utilizzato

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
    MATRIX region;

    //STRUTTURE OUTPUT MODIFICABILI
    int *QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA; //numero di risposte alle query
} params;

struct kdtree_node
{
    float medianCoordinate; //coordinata cut per il punto mediano
    float h_min, h_max;
    int indexMedianPoint;
    int numPoint;
    struct kdtree_node *left, *right;
};

/*
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
*/

float euclideanDistance(float *point, MATRIX q, int indexQ, int k)
{

    float sum = 0;
    int i;
    for (i = 0; i < k; i++)
    {
        sum += powf(point[i] - q[indexQ * k + i], 2);
    }
    return sqrtf(sum);
}

float euclideanDistanceDataset(MATRIX ds, int indexMedian, MATRIX q, int indexQ, int k)
{

    float sum = 0;
    int i;
    for (i = 0; i < k; i++)
    {
        sum += powf(ds[indexMedian * k + i] - q[indexQ * k + i], 2);
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

int countPointQuery = 0; //varibile statica per accedere alla lista dei punti vicini
int indexList = 0;


/*  Questo metodo effettua la rangeQuery sui sottoalberi figli di root, la distanza della region viene effettuata solo per le coordinate di cut del nodo padre
*   poichè le restanti coordinate non cambiano e alla fine si eseguono le chiamate ricorsive sui figli finchè non si arriva a un nodo foglia
*/
int *rangeQueryChild(MATRIX ds, KDTREE node, int n, int k, MATRIX q, int indexQ, int r, int liv, float *point, int *list)
{
    int cut = (liv - 1) % k; //il cut del nodo padre

    if (distanceChild(node, q, indexQ, k, point, cut) > r)
    {
        return list;
    }

    if (euclideanDistanceDataset(ds, node->indexMedianPoint, q, indexQ, k) <= r)
    {
        list[indexList] = node->indexMedianPoint;

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
        list[indexList] = root->indexMedianPoint;
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

// Calcola il risultato come una matrice di nQA coppie di interi
// (id_query, id_vicino)
// o in altro formato
// -------------------------------------------------

void range_query(params *input)
{
    printf("\n inizio QUERY");

    float *point = malloc(input->k * sizeof(float)); //punto che conterrà il punto più vicino alla region per il confronto
    input->QA = malloc(input->nq * sizeof(int *));
    if (point == NULL || input->QA == NULL)
    {
        printf("\nNO MEMORIA");
        exit(1);
    }
    int i, indexList;
    for (i = 0; i < input->nq; i++)
    {
        int *list = malloc(input->n * sizeof(int));
        rangeQueryRoot(input->ds, input->kdtree, input->n, input->k, input->qs, i, input->r, 0, input->region, point, list);
        // input->nQA += indexList;
        list[indexList] = -1;
        indexList = 0;
        input->QA[i] = list;

    }


}

/*
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
    // input->kdtree->point = calloc(sizeof(int), k);

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

    // if (input->h > 0)
    // {
    //     t = clock();
    //     pca(input);
    //     t = clock() - t;
    //     time = ((float)t) / CLOCKS_PER_SEC;
    //     sprintf(fname, "%s.U", input->filename);
    //     save_data(fname, input->U, input->n, input->h);
    //     sprintf(fname, "%s.V", input->filename);
    //     save_data(fname, input->V, input->k, input->h);
    // }
    // else
    //     time = -1;

    // if (!input->silent)
    //     printf("\nPCA time = %.3f secs\n", time);
    // else
    //     printf("%.3f\n", time);

    //
    // Costruzione K-d-Tree
    //

    // if (input->kdtree)
    // {
    //     t = clock();
    //     kdtree(input);
    //     t = clock() - t;
    //     time = ((float)t) / CLOCKS_PER_SEC;
    // }
    // else
    //     time = -1;
    // if (!input->silent)
    //     printf("\nIndexing time = %.3f secs\n", time);
    // else
    //     printf("%.3f\n", time);

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

    // if (input->r >= 0)
    // {
    //     if (!input->silent && input->display)
    //     {
    //         //NB: il codice non assume che QA sia ordinata per query, in caso lo sia ottimizzare il codice
    //         printf("\nQuery Answer:\n");
    //         for (i = 0; i < input->nq; i++)
    //         {
    //             printf("query %d: [ ", i);
    //             for (j = 0; j < input->nQA; j++)
    //                 if (input->QA[j * 2] == i)
    //                     printf("%d ", input->QA[j * 2 + 1]);
    //             printf("]\n");
    //         }
    //         printf("\n");
    //     }
    //     sprintf(fname, "%s.qa", input->filename);
    //     save_data(fname, input->QA, input->nQA, 2);
    // }

    if (!input->silent)
        printf("\nDone.\n");

    return 0;
}
*/
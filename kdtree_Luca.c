
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

typedef struct
{
    float *head;
    float k_min, k_max;
    struct kdtree_node *left, *right;
} kdtree_node;

#define MATRIX float *
#define KDTREE kdtree_node * // modificare con il tipo di dato utilizzato

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

kdtree_node *builtree(float *, int, int, int, int);
void order(float *, int, int, int, int);

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

// PROCEDURE ASSEMBLY
extern void prova(params *input);

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

/*
*	K-d-Tree
* 	======================
*/
//(float *ds, int liv, int k, int head, int tail)
void kdtree(params *input)
{ //provalo con valori uguali al mediano
    input->kdtree = builtree(input->ds, 0, input->k, 0, (input->k * input->n) - 1);
    //stampalbero(input->kdtree,input->k);
}

void stampalbero(kdtree_node *node, int k)
{ //da eliminare
    if (node == NULL)
    {
        printf("NULL ;\n");
        return;
    };
    for (int i = 0; i < k; i++)
    {
        printf("%f ;", node->head[i]);
    }
    printf("\n");
    stampalbero(node->left, k);
    stampalbero(node->right, k);
}

void creads(float *ds, int size)
{ //da eliminare
    for (int i = 0; i < size; i++)
    {
        ds[i] = (float)(rand() % 51) / 10;
    }
}

void stampads(float *ds, int head, int tail)
{ //da eliminare
    for (int i = head; i <= tail; i++)
        printf("_%f:", ds[i]);
    printf("\n");
}

void verificaord(float *ds, int size, int k)
{ //da eliminare
    //0 nn e' ordinato
    int tmp = 1;
    for (int i = 0; i < (size / k) - 1; i++)
        if (ds[i * k] > ds[(i + 1) * k])
        {
            tmp = 0;
            break;
        };
    printf("controllo: %i\n", tmp);
}

kdtree_node *builtree(float *ds, int liv, int k, int head, int tail)
{
    /*crea un albero con il concetto della programmazione ad oggetti, memorizzando il nodo corrente su un vettori e i riferimenti ai figli
    */
    kdtree_node *node = malloc(sizeof(kdtree_node));
    float *head_n = malloc(k * sizeof(float));
    node->head = head_n;
    if (node == NULL || head_n == NULL)
    {
        printf("Memoria insufficiente per allocare il nodo\n");
        exit(-1);
    }
    int n = ((tail - head) + 1) / k;
    if (n <= 0)
    {
        free(node);
        free(head_n);
        return NULL;
    }
    order(ds, liv % k, k, head, tail);
    node->k_min = ds[head + liv];
    node->k_max = ds[tail - k + 1 + liv];
    if (n % 2 == 0)
    {
        int s = 0;
        while (s < (n / 2) - 1 && ds[(head + ((n / 2) * k) - 2 * k * (s + 1))] >= ds[(head + ((n / 2) * k) - k * (s + 1))])
        {
            s++;
        }
        for (int i = 0; i < k; i++)
        {
            head_n[i] = ds[head + ((n / 2) * k) - k + i - s * k];
        }
        node->left = builtree(ds, (liv + 1), k, head, head + ((n / 2) * k) - k - 1 - s * k);
        node->right = builtree(ds, (liv + 1), k, head + ((n / 2) * k) - s * k, tail);
    } //if
    else
    {
        if (tail - head + 1 == k)
        {
            node->left = NULL;
            node->right = NULL;
            for (int i = 0; i < k; i++)
                head_n[i] = ds[head + i];
        }
        else
        {
            int s = 0;
            while (s < n / 2 && ds[(head + ((n / 2) * k) - k * (s + 1))] >= ds[(head + ((n / 2) * k) - (s * k))])
            {
                s++;
            }
            for (int i = 0; i < k; i++)
            {
                head_n[i] = ds[head + ((n / 2) * k) + i - s * k];
            }

            node->left = builtree(ds, (liv + 1), k, head, head + ((n / 2) * k) - 1 - s * k);
            node->right = builtree(ds, (liv + 1), k, head + ((n / 2) * k) + k - s * k, tail);
        }
    }
    return node;
}

void order(float *ds, int liv, int k, int head, int tail)
{ //implementa di k i
    int n = ((tail - head) + 1) / k;
    float *tmp = malloc(k * sizeof(float));
    for (int i = head; i < tail; i = i + k)
    {
        int idmin = i;
        float min = ds[i + liv];
        for (int j = i; j < tail; j = j + k)
        {
            if (ds[j + liv] < min)
            {
                min = ds[j + liv];
                idmin = j;
            } //if
        }     //for

        for (int g = 0; g < k; g++)
        {
            tmp[g] = ds[g + i];
            ds[g + i] = ds[g + idmin];
            ds[g + idmin] = tmp[g];
        }
    } //for
    free(tmp);
}

/*
*	Range Query Search
* 	======================
*/
void range_query(params *input)
{
    //input->qs;queryset che usa k e usa input->nq;
    //input->r;
    // for(int i=0;i<input->nq;i++){
    //rangequery(input->kdtree,);
    //  }

    // -------------------------------------------------
    // Codificare qui l'algoritmo di ricerca
    // -------------------------------------------------

    // Calcola il risultato come una matrice di nQA coppie di interi
    // (id_query, id_vicino)
    // o in altro formato
    // -------------------------------------------------
}
/*
void rangequery(kdtree_node* node,kdtree_node* nodec,float* r_query,int headq,int tailq,int r){
//kdtree_node* nodec gli passo il punto
//float* r_query,int headq,int tailq il punto di query e r
//metti i casi di uscita

}




int distance(float* r_query,int headq,int tailq,float* h){//da fare
//float* r_query,int headq,int tailq e il punto di query
//float* h rappresenta la regione indicizzata da punto
//da head e tail ti puoi prendere k
return 1;
}

int euclideandistance(float * q,float * p){
return 1;
}
*/

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

    t = clock();
    kdtree(input);
    t = clock() - t;
    time = ((float)t) / CLOCKS_PER_SEC;

    /*
    if(input->kdtree){
        t = clock();
        kdtree(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
    }else
        time = -1;
    */
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
                printf("query %d: [ ", i);
                for (j = 0; j < input->nQA; j++)
                    if (input->QA[j * 2] == i)
                        printf("%d ", input->QA[j * 2 + 1]);
                printf("]\n");
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
# K-D-Tree-and-NIPALS-Algorithm
Progetto - Architetture e programmazione dei sistemi di eleborazione

## Run
if you want to compile and run use:

  

    ./runkdtreepca32 D [-pca <h>] [-kdtree [-rq <r>]]

  

for 64 bit version use:

  

    ./runkdtreepca64 D [-pca <h>] [-kdtree [-rq <r>]]

  
  

- D: dataset name (without '.ds')

- pca <<h>h>: number of pca component

- kdtree: index the tree

- rq <<r>r>: range query with radius r

The dataset need to be bytes typed row major in this format:

```

  

numberOfColumns numberOfRows Point1FirstComponent Point1SecondComponent and so on..

  

```

  

you can find 2 datasets in the GitHub folder of the following dimension:

  

- 8000x128 (test1.ds)

- 5000x36 (test2.ds)

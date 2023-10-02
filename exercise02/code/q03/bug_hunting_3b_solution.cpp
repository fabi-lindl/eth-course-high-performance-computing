/* 
 * The following code resolves the inefficiencies of the bug hunting
 * problem of question 3b. 
 * Note: This code is not fully functional as the definiation of 
 * some of the given variables is missing. 
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

/* Variant 1. Nested parallism. Reduction of sync. barriers. */
void work(int i, int j);

void nesting(int n) 
{
    int i, j;
    #pragma omp parallel for
    for (i=0; i < n; i++) 
    {
        #pragma omp parallel for
        for (j=0; j < n; j++)
            work(i, j);
    }
}


/* Variant 2. Loop collapsing. */
void work(int i, int j);

void nesting(int n) 
{
    int i, j;
    #pragma omp parallel
    {
        // Collapse the two perfectly nested for loops. 
        #pragma omp for collapse(2)
        for (i=0; i < n; i++)
        {
            for (j=0; j < n; j++)
                work(i, j);
        }
    }
}

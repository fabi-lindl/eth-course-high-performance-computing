/* 
 * The following code resolves the inefficiencies of the bug hunting
 * problem of question 3b. 
 * Note: This code is not fully functional as the definiation of 
 * some of the given variables is missing. 
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N 1000

// Array of structs, defined elsewhere. 
extern struct data member[N];

// Returns 1 if member[i] is a "good" member, 0 otherwise, defined elsewhere. 
extern int is_good(int i);

int good_members(N);
int pos = 0;

void find_good_members() 
{
    #pragma omp parallel for
    for (int i=0; i < N; i++) 
    {
        if (is_good(i)) {
            #pragma omp critical
            {
                good_members[pos] = i;
                pos++;
            }
        }
    }
}

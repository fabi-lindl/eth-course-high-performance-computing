/* 
 * The following code resolves the bugs of the bug hunting problem 
 * of question 3a. 
 * Note: This code is not fully functional as the definiation of 
 * some of the given variables is missing. 
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void do_work(const float a, const float sum);
double new_value(int i);

void time_loop()
{
    float t = 0;
    float sum = 0; 
    int cnt = 0;

    #pragma omp parallel for num_threads(<reasonable_num>)
    for (int step=0; step < 100; step++)
    {

        #pragma omp parallel for num_threads(<reasonable_num>)
        for (int i=1; i < n; i++)
        {
            b[i - 1] = (a[i] + a[i - 1]) / 2;
            c[i - 1] += a[i];
        }

        #pragma omp parallel for num_threads(<reasonable_num>)
        for (int i=0; i < m; i++)
            z[i] = sqrt(b[i] + c[i]);

        #pragma omp critical
        {
            #pragma omp parallel num_threads(<reasonable_num>)
            #pragma omp for reduction (+:sum)
            for (int i=0; i < m; i++) 
                sum = sum + z[i];
            
            do_work(t, sum);
            
            t = new_value(cnt);
            cnt++;
        }
    }
}

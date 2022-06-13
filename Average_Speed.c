#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>



// uses the stationnary distribution to calculate the avergae speed
double Average_Speed(int nb_lanes, int Vmax) {	
   	FILE *file;
	file = fopen("Sparse/Sparse_matrix.pi","r");
	if (!file) exit(1);
	char line[256];
	int linenum=0;
	int i = -1; // i is the number of the state
	double average_speed = 0;
	int k;
	while(fgets(line, 256, file) != NULL && i < (nb_lanes * (Vmax+1)))
	{
        double value;
        linenum++;
        if(line[0] == '#') continue;

        if(sscanf(line, "%lf", &value) != 1)
        {
                printf("Syntax error, line %d\n", linenum);
                continue;
        }
        if( i != -1){
        	k = (i / nb_lanes);
			average_speed = average_speed + ( value * k );
		}
		i = i+1;
	}
	printf("Average speed = %lf\n",average_speed);
	return average_speed;
}




/* run this program using the console pauser or add your own getch, system("pause") or input loop */

int main(int argc, char *argv[]) {
	double x = Average_Speed(4,20);
	return 0;
}

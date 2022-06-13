#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define Max(x,y) (((x) >= (y)) ? (x) : (y))



/* study variables */

double Nblanes; // The number of manes on the road
double Vmax; // the max speed of each vehicule
int Nbstates; // the number of states of the markov chain

double Dab; // distance to the vehicule ahead
double Dc; // critical distance to overtake
double Va; // speed of the motorcycle
double Vb; // speed of the vehicule ahead
double Vc; // Confortable speed
double Dx; // minimal distance to vehicule bakcward
double Kx; // constant
double K1; // contant
double K2; // contant
double K3; // contant
double Sigma; // contant
double Dab2; // distance to the vehicule on the right
double Dab3; // distance to the vehicule on the left
double Dcl; // critical distance on the left
double Dcr; // critical distance on the right
double GammaL; // 
double GammaR; // 
double W; // Dab2 + Dab3
double p1; // probability of having an obstacle ahead
double p2; // probability of a slower vehicule ahead
double p3; // probability of having a confortable speed
double p4; // probability of room to maneuvre
double p5; // probability of going right being safer 
double p6; // probability of a threat
double poL; // probability of making a left maneuvre
double poR; // probability of making a right maneuvre
double pA2; // P(A2) prob of increasing speed without maneuvre.
double pA3; // P(A3) prob of decreasing speed without maneuvre.
double pA4; // P(A4) prob of keeping same speed with right maneuvre
double pA5; // P(A5) prob of increasing speed with right maneuvre
double pA6; // P(A6) prob of keeping the same speed with left maneuvre
double pA7; // P(A7) prob of increasing speed with left maneuvre


// returns the sign of the input
int sign(double Va) {
	if( Va >= 0){
		return 1;
	}else{
		return -1;
	}
}
// return a random double betweein min and max
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
// returs the probability of left maneuvring
double PoL(double Sigma, double K3, double Dab2,double W, double Dcl, double K2, double GammaL) {
	return Sigma * exp(K3 * (1 / W) * (Dab2 - Dcl) - K2 * GammaL) ;
}
// returs the probability of right maneuvring
double PoR(double Sigma, double K3, double Dab3,double W, double Dcr, double K2, double GammaR) {
	return Sigma * exp(K3 * (1 / W) * (Dab3 - Dcr) - K2 * GammaR) ;
}


// Calculates P1
double P1(double K1, double Dab, double Dc) {
	if(Dab > Dc){
		return  exp(-K1*(Dab - Dc));
	}else{
		return 1;
	}
}
// Calculates P2
double P2(double Va, double Vb,double Vmax) {
	double x;
	if(Va>Vb){
		x = 0.5 + (Va-Vb) / (2*Vmax);
		return x;
	}else{
		x = 0.5 - (Vb-Va) / (2*Vmax);
		return  x;
	}
	//return 0.5 * (sign( Va - Vb ) +1);
}
// Calculates P3 
double P3(double Va, double Vc) {
	double x;
	if(Vc>Va){
		x = 1 - (Vc-Va) / (Vc);
		return x;
	}else{
		x = 1 - (Va-Vc) / (Va);
		return  x;
	}
}
// Calculates P4
double P4(double PoL, double PoR) {
	return Max(PoL,PoR);
}
// Calculates P5
double P5(double PoL, double PoR) {
	if(PoR > PoL) {
		return 0.5 + (PoR - PoL) / 2; // si la proba de partir a droite est plus grande alors c plus safe a droite
	}else{
		return 0.5 - (PoL - PoR) / 2;
	} 
}
// Calculates P6
double P6(double Kx, double Dx) {
	return exp(-Kx*Dx); // probabilite of having a threat
}


// ----------------------------------------- Calcul des proba des actions ------------------------------------------
// P(A2) = (1-P3)(1-P1P2)
double PA2(double P1, double P2, double P3) {
	return (1-P3) * (1-P1*P2);
}
// P(A3) = P1P2(1-P4+P4P6)
double PA3(double P1, double P2, double P4, double P6) {
	return P1*P2*(1-P4+P4*P6);
}
// P(A4) = P1P2P3P4P5(1-P6)
double PA4(double P1, double P2, double P3, double P4, double P5, double P6) {
	return P1*P2*P3*P4*P5*(1-P6);
}
// P(A5) = P1P2(1-P3)P4P5(1-P6)
double PA5(double P1, double P2, double P3, double P4, double P5, double P6) {
	return P1*P2*(1-P3)*P4*P5*(1-P6);
}
// P(A6) = P1P2P3P4(1-P5)(1-P6)
double PA6(double P1, double P2, double P3, double P4, double P5, double P6) {
	return P1*P2*P3*P4*(1-P5)*(1-P6);
}
// P(A7) = P1P2(1-P3)P4(1-P5)(1-P6)
double PA7(double P1, double P2, double P3, double P4, double P5, double P6) {
	return P1*P2*(1-P3)*P4*(1-P5)*(1-P6);
}

// Loads the study parameters from a config file
void Load_parameters() {	
   	FILE *file;
	file = fopen("Configuration.in","r");
	if (!file) exit(1);
	char line[256];
	int linenum=0;
	double ret;
    char *ptr;
	while(fgets(line, 256, file) != NULL)
	{
        char param[256], value[256];

        linenum++;
        if(line[0] == '#') continue;

        if(sscanf(line, "%s %s", param, value) != 2)
        {
                printf("Syntax error, line %d\n", linenum);
                continue;
        }
		ret = strtod(value, &ptr);
		if (strcmp(param, "Dab") == 0) Dab = ret;
		if (strcmp(param, "Dc") == 0) Dc = ret;
		if (strcmp(param, "Va") == 0) Va = ret;
		if (strcmp(param, "Vb") == 0) Vb = ret;
		if (strcmp(param, "Vc") == 0) Vc = ret;
		if (strcmp(param, "Dx") == 0) Dx = ret;
		if (strcmp(param, "Kx") == 0) Kx = ret;
		if (strcmp(param, "K1") == 0) K1 = ret;
		if (strcmp(param, "K2") == 0) K2 = ret;
		if (strcmp(param, "K3") == 0) K3 = ret;
		if (strcmp(param, "Sigma") == 0) Sigma = ret;
		if (strcmp(param, "Dab2") == 0) Dab2 = ret;	
		if (strcmp(param, "Dab3") == 0) Dab3 = ret;
		if (strcmp(param, "Dcl") == 0) Dcl = ret;
		if (strcmp(param, "Dcr") == 0) Dcr = ret;
		if (strcmp(param, "GammaL") == 0) GammaL = ret;			
		if (strcmp(param, "GammaR") == 0) GammaR = ret;	
		if (strcmp(param, "Nblanes") == 0) Nblanes = ret;			
		if (strcmp(param, "Vmax") == 0) Vmax = ret;		
	}
	
	Nbstates = (int)Nblanes*(Vmax+1);
}
// the transition matrix Pxy.
double P[1000][1000];
// this will initialize the transition matrix with zeroes
void initialisation(){
	int i;
	int j;
	for (i = 0; i<1000 ; i++){
		for (j = 0; j < 1000;j++){
			P[i][j] = 0;
		}
	}
}
// this will fill the transition matrix with the right probabilities
void Create_Transition_matrix(){

 // this function will construct the transition matrix of the markov chain.
	int i,j;
	int M = Nblanes;
	for(i=0;i<Nbstates;i++){
		for(j=0;j<Nbstates;j++){
			// condition pour P(A5)
			if((i % M) != (M-1) && i + M <= Nbstates && (j == i + M + 1)){
				//printf("P(A5) in S%d-----> S%d\n",i,j);
				P[i][j] = pA5;
			}else{
				// condition pour P(A2)
				if(i + M < Nbstates && j == i + M){
					//printf("P(A2) in S%d-----> S%d\n",i,j);
					P[i][j] = pA2;
				}else{
					if((i % M != 0) && (i+M <= Nbstates) && (j == i + M - 1)){
					//	printf("P(A7) in S%d-----> S%d\n",i,j);
						P[i][j] = pA7;
					}else{
						if(i >= M && i == j + M){
								//printf("P(A3) in S%d-----> S%d\n",i,j);
								P[i][j] = pA3;
						}else{
							if(i >= M && ((i % M) != (M-1)) && j == i+1){
								//printf("P(A4) in S%d-----> S%d\n",i,j);
								P[i][j] = pA4;
							}else{
								if(i >= M+1 && (i % M != 0) && j == i-1){
									//printf("P(A6) in S%d-----> S%d\n",i,j);
									P[i][j] = pA6;
								}else{
									P[i][j] = 0;
								}
							}
						}
					}
				}
			}	
		}
	}
	
	
}


// returns the number of transititions in the markov chain based on the inputs
int Nb_transitions(){
	int i,j;
	int nb_transitions = 0;
	for(i=0;i<Nbstates;i++){
		for(j=0;j<Nbstates;j++){
			if( P[i][j] != 0) nb_transitions++;
		}
	}
	return nb_transitions;
}

// Struct of a transition
struct transition{
   int destination;
   double probability;
};

struct transition Transitions[1000];
// this will transform the transition matrix to a sparse matrix and write it to a file
void Tranform_to_sparse_matrix(){
	FILE *F;
	F = fopen("Sparse_matrix","w");
	if (!F) exit(1);
	int nb_transactions = Nb_transitions(); // calculate the number of transactions in the transition matrix.
	fprintf(F,"%d\n",nb_transactions);// we write the number of transactions to the file.
	fprintf(F,"%d\n",Nbstates);// we write the number of states to the file.
	int nb_transaction_par_ligne = 0;
	int i,j;
	int k = 0;
	int l;
	for(i=0;i<Nbstates;i++){
		for(j=0;j<Nbstates;j++){
			if(P[i][j] != 0){
				Transitions[k].destination = j;
				Transitions[k].probability = P[i][j];
				k = k+1;
				nb_transaction_par_ligne ++;
			}
		}
		fprintf(F,"%d ",nb_transaction_par_ligne);// we write the number of transactions to the file.
		for(l = 0;l<k;l++){
			if(l == k-1){
				fprintf(F,"%.3lf %d",Transitions[l].probability,Transitions[l].destination);
			}else{
				fprintf(F,"%.3lf %d ",Transitions[l].probability,Transitions[l].destination);	
			}
		}
		l = 0;
		k = 0;
		nb_transaction_par_ligne = 0;
		if(i != Nbstates - 1){
			fprintf(F,"\n");
		}
	}
	fclose(F);
}


int main(int argc, char *argv[]) {
	
	srand (time ( NULL));
	
	Load_parameters();
	
	initialisation();// make the transiton matrix full of Zeros
	
	W = Dab2 + Dab3;
			
		poL = PoL(Sigma, K3, Dab2, W, Dcl, K2, GammaL); 
		poR = PoR(Sigma, K3, Dab3, W, Dcr, K2, GammaR); 
		p1 = P1(K1, Dab, Dc);
		p2 = P2(Va, Vb, Vmax);
		p3 = P3(Va, Vc);
		p4 = P4(poL, poR);
		p5 = P5(poL, poR);
		p6 = P6(Kx, Dx);
		
		/* Calcule des proba P(A2),P(A3),..P(A7) ainsi que de POL et POR*/
		pA2 = PA2(p1, p2, p3);
		pA3 = PA3(p1, p2, p4, p6);
		pA4 = PA4(p1, p2, p3, p4, p5, p6);
		pA5 = PA5(p1, p2, p3, p4, p5, p6);
		pA6 = PA6(p1, p2, p3, p4, p5, p6);
		pA7 = PA7(p1, p2, p3, p4, p5, p6);
		printf("P1 = %.3lf | P(A2) = %.3lf | P(A3) = %.3lf | P(A4) = %.3lf | P(A5) = %.3lf | P(A6) = %.3lf | P(A7) = %.3lf\n",p1,pA2,pA3,pA4,pA5,pA6,pA7);
		//printf("Pol = %.3lf | Por = %.3lf | P1 = %.3lf | P2 = %.3lf | P3 = %.3lf | P4 = %.3lf | P5 = %.3lf | P6 = %.3lf\n",poL, poR, p1,p2,p3,p4,p5,p6);
	
	Create_Transition_matrix();
	
	Tranform_to_sparse_matrix();
	
	//FILE *file;
	//file = fopen("Variation_en_fonction_de_Dab3_Dab2.data","w");
	//if (!file) exit(1);
	/*while(Dab3 >= 1){
		W = Dab2 + Dab3;
		poL = PoL(Sigma, K3, Dab2, W, Dcl, K2, GammaL); 
		poR = PoR(Sigma, K3, Dab3, W, Dcr, K2, GammaR); 
		p1 = P1(K1, Dab, Dc);
		p2 = P2(Va, Vb, Vmax);
		p3 = P3(Va, Vc);
		p4 = P4(poL, poR);
		p5 = P5(poL, poR);
		p6 = P6(Kx, Dx);
		
		// Calcule des proba P(A2),P(A3),..P(A7) ainsi que de POL et POR
		pA2 = PA2(p1, p2, p3);
		pA3 = PA3(p1, p2, p4, p6);
		pA4 = PA4(p1, p2, p3, p4, p5, p6);
		pA5 = PA5(p1, p2, p3, p4, p5, p6);
		pA6 = PA6(p1, p2, p3, p4, p5, p6);
		pA7 = PA7(p1, p2, p3, p4, p5, p6);
		//printf("P1 = %.3lf | P(A2) = %.3lf | P(A3) = %.3lf | P(A4) = %.3lf | P(A5) = %.3lf | P(A6) = %.3lf | P(A7) = %.3lf\n",p1,pA2,pA3,pA4,pA5,pA6,pA7);
		printf("P1 = %.3lf | P(A2) = %.3lf | P(A3) = %.3lf | P(A4) = %.3lf | P(A5) = %.3lf | P(A6) = %.3lf | P(A7) = %.3lf\n",p1,pA2,pA3,pA4,pA5,pA6,pA7);
		fprintf(file,"%.1lf %.1lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",Dab3,Dab2,pA2,pA3,pA4,pA5,pA6,pA7);
		//Dab = Dab + 0.05;
		//Vc = Vc + 1;
		//Vc = Vc+1;
		Dab3 = Dab3 - 0.1;
		Dab2 = Dab2 + 0.1;
	}*/
	return 0;
}

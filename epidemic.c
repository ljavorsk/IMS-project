/*
 * Source code for ISA project.
 * file: epidemic.c
 * 
 * (C) Lukas Javorsky (xjavor20)
 * (C) Patrik Ondriga (xondri08) 
 */

#include <stdio.h>
#include <math.h>

#define BA 0;
#define TT 1;
#define TN 2;
#define NR 3;
#define ZA 4;
#define BB 5;
#define PO 6;
#define KE 7;

#define BETA 0.2
#define THETA 16
#define ALFA 1 // <0, 1>
#define GAMA 0.083
#define RO BETA / GAMA

#define DISTRICT_NUMBER 8

int m[DISTRICT_NUMBER][DISTRICT_NUMBER] = {
	{0, 43100, 12500, 27300, 13400, 15500, 11600, 10900},
	{8600, 0, 8700, 29100, 3400, 2900, 2500, 3800},
	{3900, 15400, 0, 12700, 10300, 4500, 4500, 3000},
	{3200, 16100, 11000, 0, 3600, 11400, 2000, 7300},
	{3200, 2300, 20800, 5200, 0, 9300, 10000, 6300},
	{5900, 2400, 6300, 10700, 10900, 0, 7000, 9000},
	{2500, 1200, 700, 1700, 9200, 2700, 0, 39100},
	{2800, 1400, 3400, 2400, 4600, 9800, 83400, 0}
};
int n[DISTRICT_NUMBER] = {669592, 564917, 584569, 674306, 691509, 645276, 826244, 801460};
int s[DISTRICT_NUMBER];
int I[DISTRICT_NUMBER] = {305, 315, 444, 292, 585, 274, 598, 333};
int r[DISTRICT_NUMBER] = {0, 0, 0, 0, 0, 0, 0, 0};

double getX(int district){
	return I[district] / n[district];
}

double getY(int district){
	return s[district] / n[district];
}

int sumMlj(int district){
	int result = 0;
	for(int i=0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += m[i][district];
		}
	}
	return result;
}

int sumMjl(int district){
	int result = 0;
	for(int i=0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += m[district][i];
		}
	}
	return result;
}

double sumMXBETA(int district){
	double result = 0;
	for(int i=0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += m[district][i]*getX(i)*BETA;
		}
	}
	return result;
}

double sumInfectOut(int district){
	double result = 0;
	for (int i = 0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += (m[i][district]*((I[i] - getX(i) * sumMlj(i))*BETA+sumMXBETA(i)))/(n[i]-sumMlj(i)+sumMjl(i));
		}
	}
	return result;
}

int susscpectNumber(int district){
	double nowS = s[district] - THETA * ((s[district]*I[district]*BETA )/n[district]);
	nowS -= ALFA*(1-THETA)*(((s[district]-getY(district)*sumMlj(district)) * (sumMXBETA(district) + (I[district]-getX(district)*sumMlj(district))*BETA))/(n[district]-sumMlj(district)+sumMjl(district)));
	nowS -= ALFA*(1-THETA)*getY(district)*sumInfectOut(district);
	return round(nowS);
}

int infectionNumber(int district){
	double nowI = I[district] + THETA * ((s[district]*I[district]*BETA )/n[district]) - GAMA * I[district];
	nowI += ALFA*(1-THETA)*(((s[district]-getY(district)*sumMlj(district)) * (sumMXBETA(district) + (I[district]-getX(district)*sumMlj(district))*BETA))/(n[district]-sumMlj(district)+sumMjl(district)));
	nowI += ALFA*(1-THETA)*getY(district)*sumInfectOut(district);
	return round(nowI);
} 

int main( int argc, char* argv[] )
{
	for(int i=0; i < DISTRICT_NUMBER; i++){
		s[i] = n[i] - I[i];
	}

	printf("R0: %f\n\n", RO);

	int addCases = 0;
	for(int i=0; i < 30; i++){
		int allCases = 0;
		int allSus = 0;
		int allCure = 0;
		for(int j=0; j<DISTRICT_NUMBER; j++){
			allCases += I[j];
			allSus += s[j];
			allCure += r[j];
		}
		printf("DAY %d:\t|\tnakazeni: %d   \t|\tzdravi: %d\t|\tvylieceni: %d    \t|\tpribudlo: %d     \t|\tspolu: %d\n", i, allCases, allSus, allCure, addCases, allCases+allCure+allSus);
		addCases = 0;
		for(int j=0; j<DISTRICT_NUMBER; j++){
			int newS = susscpectNumber(j);
			addCases -= I[j];
			r[j] += GAMA * I[j];
			I[j] = infectionNumber(j);
			s[j] = newS;
			addCases += I[j];
		}
	}
}

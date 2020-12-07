/*
 * Source code for IMS project.
 * file: epidemic.c
 * 
 * (C) Lukas Javorsky (xjavor20)
 * (C) Patrik Ondriga (xondri08) 
 */

#include <stdio.h>
#include <math.h>

/* Districts */
#define BA 0
#define TT 1
#define TN 2
#define NR 3
#define ZA 4
#define BB 5
#define PO 6
#define KE 7
/* --------- */

const double BETA = 0.03405;
const double THETA = 16;
const double ALFA = 0.9; // <0, 1>
const double GAMMA = 0.01775;

/* R0 = (BETA/GAMMA) 1.93 -> Reproduction number */

#define DISTRICT_NUMBER 8

/* 
	Matrix for number of people transitioning from one district to another 
	Index of district is defined above 
*/
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

/* Population of each district */
int n[DISTRICT_NUMBER] = {669592, 564917, 584569, 674306, 691509, 645276, 826244, 801460};

/* Healthy people (not infected yet) */
int s[DISTRICT_NUMBER];

/* Array of infected people in each district */
int I[DISTRICT_NUMBER] = {305, 315, 444, 292, 585, 274, 598, 333};

/* Array of recovered people in each district */
int r[DISTRICT_NUMBER] = {346, 357, 504, 332, 664, 311, 678, 379};

/**
 * Return infected ration in district.
 * @param district specific district.
 * @return double value of infected ration.
 */
double getX(int district){
	return I[district] / n[district];
}

/**
 * Return healthy people (not infected yet) ration in district.
 * @param district specific district.
 * @return double value of healthy people ration.
 */
double getY(int district){
	return s[district] / n[district];
}

/**
 * Calculate number of people, that come to this district from other.
 * @param district specific district.
 * @return number of people, that come here
 */
int sumMlj(int district){
	int result = 0;
	for(int i=0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += m[i][district];
		}
	}
	return result;
}

/**
 * Calculate number of people, that come to other districts from there.
 * @param district specific district.
 * @return number of people, that come to other districts
 */
int sumMjl(int district){
	int result = 0;
	for(int i=0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += m[district][i];
		}
	}
	return result;
}

/**
 * Helpful function for calculate sum of variables m, x and beta for all districts.
 * @param district specific district.
 * @return sum of variables m, x and beta for all districts
 */
double sumMXBETA(int district){
	double result = 0;
	for(int i=0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			result += m[district][i]*getX(i)*BETA;
		}
	}
	return result;
}

/**
 * Helpful function for calculate number of infected people, which was in other district.
 * @param district specific district.
 * @return number of infected people, which was in other district
 */
double sumInfectOut(int district){
	double result = 0;
	for (int i = 0; i < DISTRICT_NUMBER; i++){
		if(i != district){
			double var0 = I[i] - getX(i) * sumMlj(i);
			double var1 = var0 * BETA + sumMXBETA(i);
			result += m[i][district] * var1 / (n[i] - sumMlj(i) + sumMjl(i));
		}
	}
	return result;
}

/**
 * Calculate number of healthy people (not infected yet) for actual day.
 * @param district specific district.
 * @return number of healthy people (not infected yet) for actual day.
 */
int suspectNumber(int district){
	double nowS = s[district] - THETA * s[district] * I[district] * BETA / n[district];
	nowS -= ALFA*(1-THETA)*(((s[district]-getY(district)*sumMlj(district)) * 
		(sumMXBETA(district) + (I[district]-getX(district)*sumMlj(district))*BETA)) /
		(n[district]-sumMlj(district)+sumMjl(district)));
	nowS -= ALFA * (1-THETA) * getY(district) * sumInfectOut(district);
	return round(nowS);
}

/**
 * Calculate number of infected people for actual day.
 * @param district specific district.
 * @return number of infected people for actual day.
 */
int infectionNumber(int district){
	double nowI = (I[district] + THETA * s[district] * I[district] * BETA / n[district]) - GAMMA * I[district];
	nowI += ALFA*(1-THETA)*(((s[district]-getY(district)*sumMlj(district)) * 
		(sumMXBETA(district) + (I[district]-getX(district)*sumMlj(district))*BETA)) /
		(n[district]-sumMlj(district)+sumMjl(district)));
	nowI += ALFA * (1-THETA) * getY(district) * sumInfectOut(district);
	return round(nowI);
} 

int main( int argc, char* argv[] )
{
	for(int i=0; i < DISTRICT_NUMBER; i++){
		s[i] = n[i] - I[i] - r[i];
	}

	// At day 0 -> 56 cases were reported as new
	int newCases = 56;

	for(int i=0; i <= 60; i++){
		int curretlyInfected = 0;
		int currentlyHealthy = 0;
		int curedCases = 0;

		for(int j=0; j<DISTRICT_NUMBER; j++){
			curretlyInfected += I[j];
			currentlyHealthy += s[j];
			curedCases += r[j];
		}
		
		
		if(i%5 == 0){
			printf("\n--------------------------------------------------------------------------------\n");
			printf("DISTRICT\t|\tINFECTED\t|\tHEALTHY  \t|\tCURED\n");
			printf("--------------------------------------------------------------------------------\n");
			printf("Bratislava\t|\t%d       \t|\t%d    \t|\t%d\n", I[BA], s[BA], r[BA]);
			printf("Trnava   \t|\t%d       \t|\t%d    \t|\t%d\n", I[TT], s[TT], r[TT]);
			printf("Trencin  \t|\t%d       \t|\t%d    \t|\t%d\n", I[TN], s[TN], r[TN]);
			printf("Nitra    \t|\t%d       \t|\t%d    \t|\t%d\n", I[NR], s[NR], r[NR]);
			printf("Zilina   \t|\t%d       \t|\t%d    \t|\t%d\n", I[ZA], s[ZA], r[ZA]);
			printf("Banska Bystrica\t|\t%d       \t|\t%d    \t|\t%d\n", I[BB], s[BB], r[BB]);
			printf("Presov   \t|\t%d       \t|\t%d    \t|\t%d\n", I[PO], s[PO], r[PO]);
			printf("Kosice   \t|\t%d       \t|\t%d    \t|\t%d\n", I[KE], s[KE], r[KE]);
			printf("--------------------------------------------------------------------------------\n");
		}
		printf("DAY %d: | Infected: %d | Healthy: %d | Cured: %d | New Cases: %d\n", i, curretlyInfected, currentlyHealthy, curedCases, newCases);
		if(i%5 == 0){
			printf("\n");
		}
		newCases = 0;

		for(int j=0; j<DISTRICT_NUMBER; j++){
			int newS = suspectNumber(j);
			newCases -= I[j];
			r[j] += GAMMA * I[j];
			I[j] = infectionNumber(j);
			s[j] = newS;
			newCases += I[j];
		}
	}
}

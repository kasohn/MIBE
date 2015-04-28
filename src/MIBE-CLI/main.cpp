#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <random>

#pragma warning(disable:4996)

#define TO_ZERO 0.000000001
#define LN2 0.693147181

#define SIZE	2048

typedef struct _SNP{
	int id, cluster;
	char* genotype;
	double entropy, sum, maf, x2, mi;
}SNP;
typedef struct _CENTROID{
	int id, index;
	int *set;
	double entropy, ave;
}CENTROID;
typedef struct _COMBINATION
{
	int snp[5];			//maximum five-locus
	double value[2];	//0:mi, 1:chi-square
}COMBINATION;
typedef struct _BOOLEAN_EXPRESSION
{
	int *snp, *s_gate, *i_gate, *order;
	int e_gate, TP, TN, FP, FN;
	double mi;
}BOOLEAN_EXPRESSION;

//global variables
int num_case;
int num_control;
int num_sample;
int num_SNP_orin; 
int num_SNP = 0, set_size = 10, best_size = 2, from, to;
double h_type;
SNP *snp;
CENTROID centroid[5];

char f_case[200], f_control[200];
static double igf(double S, double Z)
{
	if(Z < 0.0)
	{
		return 0.0;
	}
	double Sc = (1.0 / S);
	Sc *= pow(Z, S);
	Sc *= exp(-Z);

	double Sum = 1.0;
	double Nom = 1.0;
	double Denom = 1.0;

	for(int I = 0; I < 200; I++)
	{
		Nom *= Z;
		S++;
		Denom *= S;
		Sum += (Nom / Denom);
	}

	return Sum * Sc;
}
double approx_gamma(double Z)
{
	const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
	const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

	double D = 1.0 / (10.0 * Z);
	D = 1.0 / ((12 * Z) - D);
	D = (D + Z) * RECIP_E;
	D = pow(D, Z);
	D *= sqrt(TWOPI / Z);

	return D;
} 
double chisqr(int Dof, double Cv)
{
	if(Cv < 0 || Dof < 1)
	{
		return 0.0;
	}
	double K = ((double)Dof) * 0.5;
	double X = Cv * 0.5;
	if(Dof == 2)
	{
		return exp(-1.0 * X);
	}

	double PValue = igf(K, X);
	/*if(isnan(PValue) || isinf(PValue) || PValue <= 1e-8)
	  {
	  return 1e-14;
	  } */

	PValue /= approx_gamma(K);


	return (1.0 - PValue);
}
int combin(int n, int k){
	int result=1;
	for(int i=0;i<k;i++){
		result *= (n-i);
		result /= (i+1);
	}
	return result;
}

void init_SNP(SNP *snp, int num_sample){
	snp->genotype = (char *)malloc(num_sample*sizeof(char));
	snp->entropy = 0;
	snp->sum = 0;
	snp->maf = 0;
	snp->x2 = 0;
	snp->mi = 0;
}
void init_CENTROID(CENTROID *centroid, int num_set){
	centroid->set = (int *)malloc(num_set*sizeof(int));
	centroid->entropy = 0;
	centroid->ave = 0;
}
void init_BOOLEAN_EXPRESSION(BOOLEAN_EXPRESSION *be ,int num_snp)
{
	int i;
	be->snp =(int *)malloc(num_snp*sizeof(int));
	be->s_gate =(int *)malloc(num_snp*sizeof(int));
	be->i_gate =(int *)malloc((num_snp-1)*sizeof(int));
	be->order = (int *)malloc((num_snp-1)*sizeof(int));
	for (i = 0; i < num_snp - 1; i++)
	{
		be->s_gate[i] = 0;
		be->i_gate[i] = 0;
		be->order[i] = i;
	}
	be->s_gate[num_snp - 1] = 0;
}
void make_string(char *str, BOOLEAN_EXPRESSION be, int num_snp){
	char **be_snp, s_temp[400], s_temp2[20];
	int i;
	double accuracy, sensitivity, specificity, balanced_accuracy;

	be_snp = (char **)malloc(num_snp*sizeof(char*));
	for(i=0; i<num_snp;i++)
		be_snp[i] = (char *)malloc(200*sizeof(char));


	for (int i = 0; i < num_snp; i++)
	{
		strcpy(be_snp[i], "(");
		if (0 == be.s_gate[i]) //and
		{
			//strcpy(s_temp, "a");
			s_temp[0] = (char)('a'+i);
			s_temp[1] = (char)('a'+i);
			s_temp[2] = (char)('\0');            
			strcat(be_snp[i], s_temp);
		}
		else
		{
			s_temp[0] = (char)('A'+i);
			s_temp[1] = (char)('a'+i);
			s_temp[2] = (char)('\0');
			strcat(be_snp[i], s_temp);
			strcpy(s_temp, " OR ");
			strcat(be_snp[i], s_temp);
			s_temp[0] = (char)('a'+i);
			s_temp[1] = (char)('a'+i);
			s_temp[2] = (char)('\0');    
			strcat(be_snp[i], s_temp);
		}

		strcpy(s_temp, ")");
		strcat(be_snp[i], s_temp);
	}
	for (int i = 0; i < num_snp - 1; i++)
	{
		if (0 == be.i_gate[i]) //and
		{
			strcpy(s_temp, "(");
			strcat(s_temp, be_snp[i]);
			strcat(s_temp, " AND ");
			strcat(s_temp, be_snp[be.order[i]+1]);
			strcat(s_temp, ")");
			strcpy(be_snp[be.order[i]+1], s_temp);            
		}
		else                //or
		{
			strcpy(s_temp, "(");
			strcat(s_temp, be_snp[i]);
			strcat(s_temp, " OR ");
			strcat(s_temp, be_snp[be.order[i]+1]);
			strcat(s_temp, ")");
			strcpy(be_snp[be.order[i]+1], s_temp);
		}
	}

	if (1 == be.e_gate){		
		//be_snp[num_snp - 1] += "=0";
		//strcpy(s_temp, "=0");
		strcat(be_snp[num_snp - 1], "=0");
	}
	else{
		//be_snp[num_snp - 1] += "=1";
		//strcpy(s_temp, "=1");
		strcat(be_snp[num_snp - 1], "=1");
	}

	//s_temp = "\t" + be_snp[num_snp - 1];
	strcpy(s_temp, "\t");
	strcat(s_temp, be_snp[num_snp-1]);
	accuracy = (double)(be.TP + be.TN) / (be.TP + be.TN + be.FP + be.FN);
	sensitivity = (double)(be.TP)/(be.TP+be.FN);
	specificity = (double)(be.TN)/(be.TN+be.FP);
	balanced_accuracy = (sensitivity+specificity)/2;
	//s_temp += "\t" + accuracy;
	sprintf(s_temp2, "\t%lf", accuracy);
	strcat(s_temp, s_temp2); 
	//s_temp += "\t" + sensitivity;
	sprintf(s_temp2, "\t%lf", sensitivity);
	strcat(s_temp, s_temp2);
	//s_temp += "\t" + specificity;
	sprintf(s_temp2, "\t%lf", specificity);
	strcat(s_temp, s_temp2);
	//s_temp += "\t" + balanced_accuracy;
	sprintf(s_temp2, "\t%lf", balanced_accuracy);
	strcat(s_temp, s_temp2);
	//s_temp += "\t" + be.TP;
	sprintf(s_temp2, "\t%d", be.TP);
	strcat(s_temp, s_temp2);
	//s_temp += "\t" + be.FP;
	sprintf(s_temp2, "\t%d", be.FP);
	strcat(s_temp, s_temp2);
	//s_temp += "\t" + be.FN;
	sprintf(s_temp2, "\t%d", be.FN);
	strcat(s_temp, s_temp2);
	//s_temp += "\t" + be.TN;
	sprintf(s_temp2, "\t%d",be.TN);
	strcat(s_temp, s_temp2);
	strcpy(str, s_temp);
	for(i=0; i<num_snp;i++)
		free(be_snp[i]);
	free(be_snp);
}
int and_gate(int x, int y)
{
	if (1 == x && 1 == y)
		return 1;
	else
		return 0;
}
int or_gate(int x, int y)
{
	if (1 == x || 1 == y)
		return 1;
	else
		return 0;
}
int be_run(BOOLEAN_EXPRESSION *be,int num_snp)
{
	int i;
	for (i = 0; i < num_snp; i++)       //each SNP
	{
		if (0 == be->s_gate[i])  //and
		{
			if (0 == be->snp[i])
				be->snp[i] = and_gate(0, 0);
			if (1 == be->snp[i])
				be->snp[i] = and_gate(0, 1);
			if(2 == be->snp[i])
				be->snp[i] = and_gate(1,1);
		}
		else        //or
		{
			if (0 == be->snp[i])
				be->snp[i] = or_gate(0, 0);
			if (1 == be->snp[i])
				be->snp[i] = or_gate(0, 1);
			if (2 == be->snp[i])
				be->snp[i] = or_gate(1, 1);
		}
	}
	for (i = 0; i < num_snp - 1; i++)
	{
		if (0 == be->i_gate[i]) //and
		{
			be->snp[be->order[i] + 1] = and_gate(be->snp[i], be->snp[be->order[i] + 1]);
		}
		else                //or
		{
			be->snp[be->order[i] + 1] = or_gate(be->snp[i], be->snp[be->order[i] + 1]);
		}
	}
	return be->snp[num_snp-1];
}
int next_order(BOOLEAN_EXPRESSION *be, int num_snp)
{
	int index = 0;
	while (true)
	{
		if (num_snp > index + be->order[index]+2)
		{
			be->order[index]++;
			break;
		}
		else
		{
			be->order[index] = 0;
			index++;
			if (num_snp-1 == index)
				return 1;
		}
	}
	return 0;
}
int next_s_gate(BOOLEAN_EXPRESSION *be, int num_snp)
{
	int index = 0;
	while (true)
	{
		if (0 == be->s_gate[index])
		{
			be->s_gate[index] = 1;
			break;
		}
		else
		{
			be->s_gate[index] = 0;
			index++;
			if (num_snp == index)
				return 1;
		}
	}
	return 0;
}
int next_i_gate(BOOLEAN_EXPRESSION *be, int num_snp)
{
	int index = 0;
	while (true)
	{
		if (0 == be->i_gate[index])
		{
			be->i_gate[index] = 1;
			break;
		}
		else
		{
			be->i_gate[index] = 0;
			index++;
			if (num_snp-1 == index)
				return 1;
		}
	}
	return 0;
}
int next(BOOLEAN_EXPRESSION *be,int num_snp)
{
	if (1 == next_s_gate(be, num_snp))
		if (1 == next_i_gate(be, num_snp))
			if (1 == next_order(be, num_snp))
				return 1;
	return 0;
}
void be_copy(BOOLEAN_EXPRESSION *to, BOOLEAN_EXPRESSION *from, int k_num){
	for (int i = 0; i < k_num; i++)
	{
		to->s_gate[i] = from->s_gate[i];
	}
	for (int i = 0; i < k_num - 1; i++)
	{
		to->i_gate[i] = from->i_gate[i];
		to->order[i] = from->order[i];
	}
	to->e_gate = from->e_gate;
	to->TP = from->TP; to->TN = from->TN; to->FP = from->FP; to->FN = from->FN;
	to->mi = from->mi;
}
void be_free(BOOLEAN_EXPRESSION *be){
	free(be->i_gate);
	free(be->order);
	free(be->s_gate);
	free(be->snp);
}
/*const bool cmp(const SNP& a, const SNP& b)
  {
  return a.mi > b.mi;
  }

  void sort()
  {
  sort(snp, snp+num_data, cmp);
  }*/

void loadData()
{
	FILE *ifp, *ifp_case, *ifp_control;
	char s_temp[10000];
	int random_seed;
	int i, j;
	int cnt = 0, cnt_case = 0, cnt_control = 0;    
	double maf, missing;

	printf("Data loading\n");

	if(NULL == (ifp = fopen("option.txt", "r"))){
		printf("No option.txt file error");
		exit(1);
	}
	fscanf(ifp, "%s", s_temp);	//case_file:
	fscanf(ifp, "%s", f_case);
	fscanf(ifp, "%s", s_temp);	//control_file:
	fscanf(ifp, "%s", f_control);
	fscanf(ifp, "%s", s_temp);	//random_seed:
	fscanf(ifp, "%d", &random_seed);
	fscanf(ifp, "%s", s_temp);	//candidates:
	fscanf(ifp, "%d", &set_size);
	fscanf(ifp, "%s", s_temp);	//order:
	fscanf(ifp, "%d", &from);
	fscanf(ifp, "%c", &s_temp[0]);
	fscanf(ifp, "%d", &to);
	fscanf(ifp, "%s", s_temp);	//boolean_expression
	fscanf(ifp, "%d", &best_size);

	srand(random_seed);
	fclose(ifp);	

	if(NULL == (ifp_case = fopen(f_case, "r"))){
		printf("No %s file error", f_case);
		exit(1);
	}
	if(NULL == (ifp_control = fopen(f_control, "r"))){
		printf("No %s file error", f_control);
		exit(1);
	}

	i=0;
	while(1)        
	{
		char *data_case, *data_control;
		int table[] = { 0, 0, 0 };
		double p_index[] = { 0, 0, 0 };
		double p_index1[2][3] = {0};
		if (NULL == fgets(s_temp, 10000, ifp_case))
			break;
		num_case = (int)strlen(s_temp)-1;		// -1 because of \n 
		data_case = (char*)malloc(num_case * sizeof(char));

		for (j = 0; j < num_case; j++)
			data_case[j] = (char)s_temp[j];
		fgets(s_temp, 10000, ifp_control);
		num_control = (int)strlen(s_temp)-1;		// -1 because of \n 
		data_control = (char *)malloc(num_control * sizeof(char));

		num_sample = num_case + num_control;

		h_type = -(double)num_case / num_sample *log((double)num_case / num_sample) / log(2.0)
			- (double)num_control / num_sample * log((double)num_control / num_sample) / log(2.0);

		for (j = 0; j < num_control; j++)
			data_control[j] = (char)s_temp[j];

		//cal_entropy        
		cnt=0;
		for (j = 0; j < num_case; j++)
		{
			if ('0' <= data_case[j] && data_case[j] <= '2')
			{
				table[data_case[j] - '0']++;
				p_index1[0][data_case[j] - '0']++;
				cnt++;
				cnt_case++;
			}
		}
		for (j = 0; j < num_control; j++)
		{
			if ('0' <= data_control[j] && data_control[j] <= '2')
			{
				table[data_control[j] - '0']++;
				p_index1[1][ data_control[j] - '0']++;
				cnt++;
				cnt_control++;
			}
		}
		for (j = 0; j < 3; j++)
		{
			if (table[j] > 0)
				p_index[j] = (double)table[j] / cnt;
			else
				p_index[j] = TO_ZERO;

			if (p_index1[0][j] > 0)
				p_index1[0][j] /= cnt;
			else
				p_index1[0][j] = TO_ZERO;

			if (p_index1[1][j] > 0)
				p_index1[1][j] /= cnt;
			else
				p_index1[1][j] = TO_ZERO;
		}
		if (p_index[0] > p_index[2])
			maf = 2 * p_index[2] + p_index[1];
		else
			maf = 2 * p_index[0] + p_index[1];
		maf /= 2;
		missing = (double)(num_sample - cnt) / num_sample;
		if (maf > 0.05 && missing < 0.05)
		{
			init_SNP(&snp[num_SNP], num_sample);
			snp[num_SNP].entropy = 0;
			snp[num_SNP].mi = 0;
			snp[num_SNP].id = i+1;             

			for (j = 0; j < num_case; j++)
				snp[num_SNP].genotype[j] = data_case[j];
			for (j = 0; j < num_control; j++)
				snp[num_SNP].genotype[num_case + j] = data_control[j];
			for (j = 0; j < 3; j++)
			{
				snp[num_SNP].entropy += (-(p_index[j]) * (log(p_index[j]) / log(2.0)));
				snp[num_SNP].mi += (p_index1[0][j] * (log(p_index1[0][j]) / log(2.0)));
				snp[num_SNP].mi += (p_index1[1][j] * (log(p_index1[1][j]) / log(2.0)));
			}
			snp[num_SNP].mi += (snp[num_SNP].entropy + h_type);
			snp[num_SNP].maf = maf;
			num_SNP++;

			if (0 == num_SNP % SIZE){
				printf("%d SNPs are loaded\n", num_SNP);
				snp = (SNP *)realloc(snp, (num_SNP + SIZE)*sizeof(SNP));
			}
		}


		i++;
	}

	fclose(ifp_case);
	fclose(ifp_control);

	printf("%d SNPs are remianed\n", num_SNP);    
}
int cmp(const void* a, const void * b)
{
	if(((SNP *)a)->mi < ((SNP *)b)->mi)
		return 1;
	else if (((SNP *)a)->mi > ((SNP *)b)->mi)
		return -1;
	else
		return 0;
}
int combiCompare(const void * x, const void * y)
{
	if(((COMBINATION *)x)->value[0] < ((COMBINATION *)y)->value[0])
		return 1;
	else if(((COMBINATION *)x)->value[0] < ((COMBINATION *)y)->value[0])
		return -1;
	else
		return 0;
}
double distance_2(int a, int b)
{	
	double h_xy = 0, h_xyt=0, h_xt=0, h_yt=0;
	double table[3][3];
	double t_case[3][3];
	double t_control[3][3];
	double px_case[3], px_control[3], py_case[3], py_control[3];
	double P_AB, P_A=0, P_B=0, P_a=0, P_b=0;
	double t_ratio;
	int cnt = 0, i, j;

	for(i=0; i<3;i++){
		for(j=0;j<3;j++){
			table[i][j] = 0;
			t_case[i][j] = 0;
			t_control[i][j] = 0;
		}
		px_case[i]=0;
		px_control[i]=0;
		py_case[i]=0;
		py_control[i]=0;
	}

	for(i=0; i<num_sample; i++)
	{
		if(snp[a].genotype[i]-'0' >= 0 && snp[b].genotype[i] -'0' >=0)
		{
			table[snp[a].genotype[i]-'0'][snp[b].genotype[i]-'0']++;

			cnt++;
			if(i<num_case)
				t_case[snp[a].genotype[i]-'0'][snp[b].genotype[i]-'0']++;
			else
				t_control[snp[a].genotype[i]-'0'][snp[b].genotype[i]-'0']++;
		}
	}

	t_ratio = (double)num_case/num_control;
	for(int i=0; i<3; i++)
		for(int j=0;j<3; j++)
		{
			table[i][j] /= cnt;
			if(table[i][j] > TO_ZERO)
				h_xy += -table[i][j] * log(table[i][j])/log(2.0);

			t_case[i][j] /= cnt;
			t_control[i][j] /= cnt;
			if(t_case[i][j] > TO_ZERO)
				h_xyt += -t_case[i][j] * log(t_case[i][j])/log(2.0);
			if(t_control[i][j] > TO_ZERO)
				h_xyt += -t_control[i][j] * log(t_control[i][j])/log(2.0);						
		}

	P_AB = table[0][0] + (table[0][1] + table[1][0])/2 + table[1][1]/4;


	for(i=0;i<3; i++){
		P_A += table[0][i] + table[1][i]/2;
		P_a += table[2][i] + table[1][i]/2;
		P_B += table[i][0] + table[i][1]/2;
		P_b += table[i][2] + table[i][1]/2;
		for(j=0;j<3;j++)
		{
			px_case[i] += t_case[i][j];
			px_control[i] += t_control[i][j];
			py_case[i] += t_case[j][i];
			py_control[i] += t_control[j][i];
		}
		if(px_case[i] > TO_ZERO)
			h_xt += -px_case[i]*log(px_case[i])/log(2.0);
		if(px_control[i] > TO_ZERO)
			h_xt += -px_control[i]*log(px_control[i])/log(2.0);
		if(py_case[i] > TO_ZERO)
			h_yt += -py_case[i]*log(py_case[i])/log(2.0);
		if(py_control[i] > TO_ZERO)
			h_yt += -py_control[i]*log(py_control[i])/log(2.0);
	}	        

	return (-h_xyt+h_type+h_xy);
}
double distance_3(int a, int b, int c)
{
	double h_xy = 0, h_xyt = 0;
	double table[3][3][3];
	double t_case[3][3][3];
	double t_control[3][3][3];
	double pa_case[3], pa_control[3], pb_case[3], pb_control[3], pc_case[3], pc_control[3];           
	double t_ratio;
	int cnt = 0, cnt_case=0, cnt_control = 0;

	int i,j,k;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				table[i][j][k] = 0;
				t_case[i][j][k] = 0;
				t_control[i][j][k] = 0;
			}
		}
		pa_case[i] = 0;
		pa_control[i] = 0;
		pb_case[i] = 0;
		pb_control[i] = 0;
		pc_case[i] = 0;
		pc_control[i] = 0;
	}

	for (int i = 0; i < num_sample; i++)
	{
		if (snp[a].genotype[i] - '0' >= 0 && snp[b].genotype[i] - '0' >= 0 && snp[c].genotype[i] - '0' >= 0)
		{
			table[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0']++;
			cnt++;
			if (i < num_case)
			{
				t_case[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0']++;
				cnt_case++;
			}
			else
			{
				t_control[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0']++;
				cnt_control++;
			}
		}
	}

	t_ratio = (double)num_case / num_control;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			for(k=0; k<3; k++)
			{
				table[i][j][k] /= cnt;
				if (table[i][j][k] > TO_ZERO)
					h_xy += -table[i][j][k] * log(table[i][j][k]) / log(2.0);

				t_case[i][j][k] /= cnt;
				t_control[i][j][k] /= cnt;
				if (t_case[i][j][k] > TO_ZERO)
					h_xyt += -t_case[i][j][k] * log(t_case[i][j][k]) / log(2.0);
				if (t_control[i][j][k] > TO_ZERO)
					h_xyt += -t_control[i][j][k] * log(t_control[i][j][k]) / log(2.0);
			}
		}

	return (-h_xyt + h_type + h_xy);
}
double distance_4(int a, int b, int c, int d)
{
	double h_xy = 0, h_xyt = 0;
	double table[3][3][3][3];
	double t_case[3][3][3][3];
	double t_control[3][3][3][3];
	double pa_case[3], pa_control[3], 
	       pb_case[3], pb_control[3],
	       pc_case[3], pc_control[3],
	       pd_case[3], pd_control[3];

	double t_ratio;
	int cnt = 0, cnt_case = 0, cnt_control = 0;

	int i, j, k, l;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				for (l = 0; l < 3; l++)
				{
					table[i][j][k][l] = 0;
					t_case[i][j][k][l] = 0;
					t_control[i][j][k][l] = 0;
				}

		pa_case[i] = 0;
		pa_control[i] = 0;
		pb_case[i] = 0;
		pb_control[i] = 0;
		pc_case[i] = 0;
		pc_control[i] = 0;
		pd_case[i] = 0;
		pd_control[i] = 0;
	}

	for (i = 0; i < num_sample; i++)
	{
		if (snp[a].genotype[i] - '0' >= 0 && snp[b].genotype[i] - '0' >= 0 && snp[c].genotype[i] - '0' >= 0 && snp[d].genotype[i] - '0' >= 0)
		{
			table[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0'][snp[d].genotype[i] - '0']++;
			cnt++;
			if (i < num_case)
			{
				t_case[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0'][snp[d].genotype[i] - '0']++;
				cnt_case++;
			}
			else
			{
				t_control[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0'][snp[d].genotype[i] - '0']++;
				cnt_control++;
			}
		}
	}

	t_ratio = (double)num_case / num_control;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				for (l = 0; l < 3; l++)
				{
					table[i][j][k][l] /= cnt;
					if (table[i][j][k][l] > TO_ZERO)
						h_xy += -table[i][j][k][l] * log(table[i][j][k][l]) / log(2.0);

					t_case[i][j][k][l] /= cnt;
					t_control[i][j][k][l] /= cnt;
					if (t_case[i][j][k][l] > TO_ZERO)
						h_xyt += -t_case[i][j][k][l] * log(t_case[i][j][k][l]) / log(2.0);
					if (t_control[i][j][k][l] > TO_ZERO)
						h_xyt += -t_control[i][j][k][l] * log(t_control[i][j][k][l]) / log(2.0);
				}                        

	return (-h_xyt + h_type + h_xy);
}
double distance_5(int a, int b, int c, int d, int e){
	double h_xy = 0, h_xyt = 0;
	double table[3][3][3][3][3];
	double t_case[3][3][3][3][3];
	double t_control[3][3][3][3][3];
	double pa_case[3], pa_control[3],
	       pb_case[3], pb_control[3],
	       pc_case[3], pc_control[3],
	       pd_case[3], pd_control[3],
	       pe_case[3], pe_control[3];

	double t_ratio;
	int cnt = 0, cnt_case = 0, cnt_control = 0;

	int i, j, k, l, m;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				for (l = 0; l < 3; l++)
					for(m=0; m<3; m++)
					{
						table[i][j][k][l][m] = 0;
						t_case[i][j][k][l][m] = 0;
						t_control[i][j][k][l][m] = 0;
					}

		pa_case[i] = 0;                pa_control[i] = 0;
		pb_case[i] = 0;                pb_control[i] = 0;
		pc_case[i] = 0;                pc_control[i] = 0;
		pd_case[i] = 0;                pd_control[i] = 0;
		pe_case[i] = 0;                pe_control[i] = 0;
	}

	for (i = 0; i < num_sample; i++)
	{
		if (snp[a].genotype[i] - '0' >= 0 && snp[b].genotype[i] - '0' >= 0 && snp[c].genotype[i] - '0' >= 0 && snp[d].genotype[i] - '0' >= 0 && snp[e].genotype[i] - '0' >= 0)
		{
			table[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0'][snp[d].genotype[i] - '0'][snp[e].genotype[i] - '0']++;
			cnt++;
			if (i < num_case)
			{
				t_case[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0'][snp[d].genotype[i] - '0'][snp[e].genotype[i] - '0']++;
				cnt_case++;
			}
			else
			{
				t_control[snp[a].genotype[i] - '0'][snp[b].genotype[i] - '0'][snp[c].genotype[i] - '0'][snp[d].genotype[i] - '0'][snp[e].genotype[i] - '0']++;
				cnt_control++;
			}
		}
	}

	t_ratio = (double)num_case / num_control;
	for ( i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				for (l = 0; l < 3; l++)
					for(m=0; m<3; m++)
					{
						table[i][j][k][l][m] /= cnt;
						if (table[i][j][k][l][m] > TO_ZERO)
							h_xy += -table[i][j][k][l][m] * log(table[i][j][k][l][m]) / log(2.0);

						t_case[i][j][k][l][m] /= cnt;
						t_control[i][j][k][l][m] /= cnt;
						if (t_case[i][j][k][l][m] > TO_ZERO)
							h_xyt += -t_case[i][j][k][l][m] * log(t_case[i][j][k][l][m]) / log(2.0);
						if (t_control[i][j][k][l][m] > TO_ZERO)
							h_xyt += -t_control[i][j][k][l][m] * log(t_control[i][j][k][l][m]) / log(2.0);
					}

	return (-h_xyt + h_type + h_xy);
}
void k_means(int k_num, int iteration)
{
	int cnt, changed, min_index=0;
	double temp_distance=0, min_distance;
	int i,j, k, t, m;

	printf("%d-locus clustering\n", k_num);        

	//Initialization
	for(i=0; i<k_num; i++){
		init_CENTROID(&centroid[i], set_size);
		centroid[i].index = rand()%num_SNP;
		centroid[i].entropy = snp[centroid[i].index].entropy;
		centroid[i].ave = 0;
	}

	changed = 1;
	for(t=0; t<iteration && changed >0 ; t++)
	{        
		changed = 0;
		//Finding the closest centroid
		for(k=0;k<k_num; k++)
			centroid[k].ave = 0;
		for(j=0; j<num_SNP; j++)
		{
			min_index = 0;	min_distance = 10000;
			for(k=0; k<k_num; k++)
			{
				temp_distance = distance_2(centroid[k].index, j);

				if(min_distance > temp_distance)
				{
					min_index = k;
					min_distance = temp_distance;
				}

			}
			if(snp[j].cluster != min_index)
			{
				changed++;
				snp[j].cluster = min_index;				
			}
			centroid[min_index].ave += min_distance;
		}

		//reclaculating centroids		
		for(m=0; m<k_num; m++)
		{
			double **genotype;
			genotype = (double **)malloc(num_sample*sizeof(double *));
			for(j=0; j<num_sample; j++)
			{
				genotype[j] = (double *)malloc(3*sizeof(double));
				genotype[j][0] = genotype[j][1] = genotype[j][2] = 0;
			}

			cnt = 0;
			for(j=0; j<num_SNP; j++)
			{
				if(m == snp[j].cluster)
				{	
					for(int k=0; k<num_sample; k++)
					{
						if(snp[j].genotype[k]-'0' >= 0)
							genotype[k][snp[j].genotype[k]-'0']++;
					}
					cnt++;
				}
			}
			for(j=0; j<num_sample; j++)
				for(k=0; k<3; k++)
				{
					genotype[j][k] /= cnt;
					if(genotype[j][k] < TO_ZERO)
						genotype[j][k] = TO_ZERO;
				}

			//대표 sNP 선정
			cnt = 0; centroid[m].entropy = snp[centroid[m].index].entropy;
			for(j=0; j<num_SNP; j++)
			{
				if(m == snp[j].cluster)
				{
					if(2 == k_num)
						snp[j].sum = distance_2(centroid[m].index, j);

					if(0 == cnt || snp[j].sum < snp[min_index].sum)
						min_index = j;

					cnt++;
				}
			}			       

			centroid[m].index = min_index;
			centroid[m].ave /= cnt;
			centroid[m].entropy = snp[min_index].entropy;
			centroid[m].id = snp[min_index].id;
		}
	}
}
void check_set(int k_num){
	int i, j;

	printf("%d-locus candidate selection\n", k_num);

	//sum of distance between SNP and other centroids
	for(i=0; i<num_SNP; i++)
	{
		snp[i].sum = 0;
		for(j=0; j<k_num; j++)
			if(j != snp[i].cluster)
				snp[i].sum += distance_2(centroid[j].index, i);
	}

	for(i=0; i<k_num ; i++)
		for(j=0; j<set_size; j++)
			centroid[i].set[j] = -1;

	//대표 SNP 집합 선택하기
	for(int i=0; i<k_num ;i++){
		for(int j=0; j<num_SNP; j++){
			if(snp[j].cluster == i){
				for(int k=0; k<set_size; k++){
					if(0 > centroid[i].set[k]){	//비어있을 경우
						centroid[i].set[k] = j;
						break;
					}
					else if(snp[centroid[i].set[k]].sum <= snp[j].sum){
						for(int l=set_size-1; l>k; l--)
							centroid[i].set[l] = centroid[i].set[l-1];	//한칸씩 뒤로 밀기
						centroid[i].set[k] = j;
						break;
					}
				}
			}
		}
	}	        
}
void single_locus()
{
	FILE *ofp;
	int i;
	char str_file[200];

	printf("Single-locus\n");
	qsort(snp, num_SNP, sizeof(SNP), cmp);

	strcpy(str_file, f_case);
	strcat(str_file, "_single.txt");
	if(NULL == (ofp = fopen(str_file, "w"))){
		printf("Cannot make %s\n", str_file);
		exit(1);
	}

	//var d = Gamma.WithShapeScale(1, 1.0/(num_sample*LN2));   //alpha, beta(1/(n*ln2))

	fprintf(ofp, "ID\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tmaf\n");


	for(i=0; i<num_SNP; i++){
		fprintf(ofp,"%d\t%lf\t%lf\t", snp[i].id, snp[i].mi, chisqr(2,2 *num_sample* LN2 *snp[i].mi)*combin(num_SNP,1));
		fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 *num_sample* LN2 *snp[i].mi,2, num_SNP, 1);
		fprintf(ofp, "\t%lf\n", snp[i].maf);
	}


	fclose(ofp);
	printf("Single-locus end\n");
}
void two_locus()
{
	int k_num = 2;
	printf("########### two-locus ###########\n");

	k_means(k_num, set_size);
	check_set(k_num);


	int *temp_set = (int *)malloc(set_size * k_num*sizeof(int));
	int temp_cnt = 0;            

	for (int i = 0; i < set_size * k_num; i++)
	{
		temp_set[i] = centroid[i / set_size].set[i % set_size];
		if (temp_set[i] >= 0)
			temp_cnt++;
	}

	int *all_combi = (int*)malloc(2*temp_cnt * (temp_cnt-1)*sizeof(int));
	double *all_value = (double*)malloc(temp_cnt * (temp_cnt-1)*sizeof(double));

	printf("two-locus best_search\n");

	int cnt = 0, total = temp_cnt * (temp_cnt-1)/2;   

	for (int i = 0; i < set_size * k_num - 1; i++)
	{        
		for (int j = i + 1; j < set_size * k_num; j++)
		{                      
			if (0 <= temp_set[i] && 0 <= temp_set[j])
			{
				all_combi[2*cnt] = temp_set[i];
				all_combi[2*cnt+1] = temp_set[j];
				all_value[cnt] = distance_2(temp_set[i], temp_set[j]);
				cnt++; 
			}
		}
	}              

	printf("two-locus sorting\n");
	//sorting////////////
	for(int i=0; i<cnt-1; i++)
		for (int j = i + 1; j < cnt; j++)
		{
			int i_temp;
			double d_temp;
			if (all_value[i] < all_value[j])
			{
				i_temp = all_combi[2*i]; all_combi[2*i] = all_combi[2*j]; all_combi[2*j] = i_temp;
				i_temp = all_combi[2*i+1]; all_combi[2*i+1] = all_combi[2*j+1]; all_combi[2*j+1] = i_temp;
				d_temp = all_value[i]; all_value[i] = all_value[j]; all_value[j] = d_temp;
			}
		}

	FILE *ofp;
	char str_file[200];
	strcpy(str_file, f_case);
	strcat(str_file, "_two_all.txt");

	ofp = fopen(str_file, "w");
	fprintf(ofp,"rank\tSNP1\tSNP2\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tboolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");


	for (int i = 0; i < cnt; i++)
	{
		fprintf(ofp, "%d\t%d\t%d\t%lf\t", i+1, snp[all_combi[2*i]].id, snp[all_combi[2*i+1]].id, all_value[i]);
		if (0 < all_value[i])
			fprintf(ofp, "%lf", chisqr(8, 2 * num_sample * LN2 * all_value[i])*combin(num_SNP,2));
		else
			fprintf(ofp, "1.0");
		fprintf(ofp, "\t");
		if (0 < all_value[i])
			fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_value[i],8, num_SNP, 2);
		//fprintf(ofp, "%lf", 2 * num_sample * LN2 * all_value[i]);
		else
			fprintf(ofp, "1.0");
		if (i < best_size)
		{
			BOOLEAN_EXPRESSION best;
			BOOLEAN_EXPRESSION be;
			best.mi = 0;
			init_BOOLEAN_EXPRESSION(&best, k_num);
			init_BOOLEAN_EXPRESSION(&be, k_num);

			int t_cnt = 0;
			while (true)        //test all boolean expression
			{
				t_cnt++;
				int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;

				for (int k = 0; k < num_sample; k++)
				{
					//label12.Text = "for test";
					for (int j = 0; j < k_num; j++)
						be.snp[j] = snp[all_combi[2*i+j]].genotype[k] - '0';
					value = be_run(&be, k_num);

					if (1 == value)
					{
						if (k < num_case)
							case_1++;
						else
							control_1++;
					}
					else
					{
						if (k < num_case)
							case_0++;
						else
							control_0++;
					}                    
				}
				if (case_0 * control_1 > case_1 * control_0)
				{
					be.e_gate = 1;
					be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
				}
				else
				{
					be.e_gate = 0;
					be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
				}
				be.mi = 0;
				be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
				be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
				be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
				be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
				be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
				be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
				//s_temp = be.mi.ToString();
				//outStream1.WriteLine(s_temp);
				if (be.mi > best.mi)
					be_copy(&best, &be, k_num);
				if (1 == next(&be, k_num))
					break;
			}

			char s_temp[200];
			make_string(s_temp, best, k_num);
			fprintf(ofp, "%s\n", s_temp); 
			be_free(&best);
			be_free(&be);
		}
		else
			fprintf(ofp, "\n"); 
	}
	fclose(ofp);
	strcpy(str_file, f_case);
	strcat(str_file, "_two_deduplicated.txt");        
	ofp = fopen(str_file, "w");

	fprintf(ofp, "rank\tSNP1\tSNP2\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_epression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

	for (int i = 0; i < cnt; i++)
	{
		int d_cnt = 0;
		for(int j=0; j<set_size*k_num;j++)
			for(int k =0; k<k_num; k++)
				if(temp_set[j] == all_combi[2*i+k])
					d_cnt++;


		if (d_cnt == k_num)
		{
			BOOLEAN_EXPRESSION best;
			BOOLEAN_EXPRESSION be;            
			init_BOOLEAN_EXPRESSION(&best, k_num);
			init_BOOLEAN_EXPRESSION(&be, k_num);
			best.mi = 0;

			int t_cnt = 0;
			while (true)        //test all boolean expression
			{
				t_cnt++;
				int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


				for (int k = 0; k < num_sample; k++)
				{
					for (int j = 0; j < k_num; j++)
						be.snp[j] = snp[all_combi[2*i+j]].genotype[k] - '0';
					value = be_run(&be, k_num);
					if (1 == value)
					{
						if (k < num_case)
							case_1++;
						else
							control_1++;
					}
					else
					{
						if (k < num_case)
							case_0++;
						else
							control_0++;
					}
					if (case_0 * control_1 > case_1 * control_0)
					{
						be.e_gate = 1;
						be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
					}
					else
					{
						be.e_gate = 0;
						be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
					}
				}
				be.mi = 0;
				be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
				be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
				be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
				be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
				be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
				be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
				if(be.mi > best.mi)
					be_copy(&best, &be, k_num);
				if (1 == next(&be, k_num))
					break;
			}	


			fprintf(ofp, "%d\t%d\t%d\t%lf\t", i+1, snp[all_combi[2*i]].id, snp[all_combi[2*i+1]].id, all_value[i]);
			if (0 < all_value[i])
				fprintf(ofp, "%lf", chisqr(8, 2 * num_sample * LN2 * all_value[i]), combin(num_SNP, 2));
			else
				fprintf(ofp, "1.0");
			fprintf(ofp, "\t");
			if (0 < all_value[i])
				fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_value[i],8, num_SNP, 2);
			//fprintf(ofp, "%lf", 2 * num_sample * LN2 * all_value[i]);
			else
				fprintf(ofp, "1.0");

			char s_temp[200];
			make_string(s_temp, best, k_num);
			fprintf(ofp, "%s\n", s_temp); 

			for (int j = 0; j < set_size * k_num; j++)
				for (int k = 0; k < k_num; k++)
					if (temp_set[j] == all_combi[2*i+k])
						temp_set[j] = -1;

			be_free(&best);
			be_free(&be);
		}
	}
	fclose(ofp);
	printf("two-locus end\n");
}
void three_locus()
{
	int k_num = 3, opt = 0;

	printf("##############   three-locus   ############\n");
	k_means(k_num, set_size);
	check_set(k_num);


	int *temp_set =  (int*)malloc(set_size * k_num*sizeof(int));
	int temp_cnt = 0;

	for (int i = 0; i < set_size * k_num; i++)
	{
		temp_set[i] = centroid[i / set_size].set[i % set_size];
		if (temp_set[i] >= 0)
			temp_cnt++;
	}

	int *all_combi = (int*)malloc(3 * temp_cnt * (temp_cnt - 1) * (temp_cnt - 2)  / 2/3 * sizeof(int));
	double* all_value = (double *)malloc((temp_cnt * (temp_cnt - 1) * (temp_cnt - 2)) / 2 / 3 * sizeof(double));

	int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2)) / 2 / 3;
	for (int i = 0; i < set_size * k_num - 2; i++)
	{       
		for (int i1 = i + 1; i1 < set_size * k_num-1; i1++)
			for(int i2 = i1+1;i2<set_size*k_num;i2++)
				if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2])
				{
					all_combi[3*cnt] = temp_set[i];
					all_combi[3*cnt+1] = temp_set[i1];
					all_combi[3*cnt+2] = temp_set[i2];
					all_value[cnt] = distance_3(temp_set[i], temp_set[i1], temp_set[i2]);                    
					cnt++;
				}
	}
	printf("three-locus sorting\n");
	//sorting////////////
	for (int i = 0; i < cnt - 1; i++)
		for (int j = i + 1; j < cnt; j++)
		{
			int i_temp;
			double d_temp;
			if (all_value[i] < all_value[j])
			{
				for (int k = 0; k < k_num; k++)
				{ i_temp = all_combi[i*3+ k]; all_combi[i*3+ k] = all_combi[3*j+k]; all_combi[3*j+k] = i_temp; }

				d_temp = all_value[i]; all_value[i] = all_value[j]; all_value[j] = d_temp;
			}

		}

	FILE *ofp;
	char str_file[200];
	strcpy(str_file, f_case);
	strcat(str_file, "_three_all.txt");

	ofp = fopen(str_file, "w");

	fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

	for (int i = 0; i < cnt; i++)
	{
		fprintf(ofp, "%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[3*i]].id, snp[all_combi[3*i+1]].id,
				snp[all_combi[3*i+ 2]].id, all_value[i]);
		if (0 < all_value[i])
			fprintf(ofp, "%lf", chisqr(26, 2 * num_sample * LN2 * all_value[i])* combin(num_SNP, 3));
		else
			fprintf(ofp, "1.0");
		fprintf(ofp, "\t");
		if (0 < all_value[i])
			fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_value[i],26, num_SNP, 3);
		else
			fprintf(ofp, "1.0");

		if (i < best_size)
		{
			BOOLEAN_EXPRESSION best;
			BOOLEAN_EXPRESSION be;

			init_BOOLEAN_EXPRESSION(&best, k_num);
			init_BOOLEAN_EXPRESSION(&be, k_num);

			best.mi = 0;

			int t_cnt = 0;
			while (true)        //test all boolean expression
			{
				t_cnt++;
				int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;

				for (int k = 0; k < num_sample; k++)
				{
					for (int j = 0; j < k_num; j++)
						be.snp[j] = snp[all_combi[3*i+j]].genotype[k] - '0';
					value = be_run(&be, k_num);
					if (1 == value)
					{
						if (k < num_case)
							case_1++;
						else
							control_1++;
					}
					else
					{
						if (k < num_case)
							case_0++;
						else
							control_0++;
					}
					if (case_0 * control_1 > case_1 * control_0)
					{
						be.e_gate = 1;
						be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
					}
					else
					{
						be.e_gate = 0;
						be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
					}
				}
				be.mi = 0;
				be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
				be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
				be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
				be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
				be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
				be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
				//s_temp = be.mi.ToString();
				//outStream1.WriteLine(s_temp);
				if (be.mi > best.mi)
					be_copy(&best, &be, k_num);
				if (1 ==next(&be, k_num))
					break;

			}		

			char s_temp[200];
			make_string(s_temp, best, k_num);
			fprintf(ofp, "%s\n", s_temp);
			be_free(&best);
			be_free(&be);
		}
		else
			fprintf(ofp,"\n");
	}
	fclose(ofp);

	strcpy(str_file, f_case);
	strcat(str_file, "_three_deduplicated.txt");

	ofp = fopen(str_file, "w");
	fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

	for (int i = 0; i < cnt; i++)
	{
		int d_cnt = 0;
		for(int j=0; j<set_size*k_num;j++)
			for(int k =0; k<k_num; k++)
				if(temp_set[j] == all_combi[2*i+k])
					d_cnt++;


		if (d_cnt == k_num)
		{
			BOOLEAN_EXPRESSION best;
			BOOLEAN_EXPRESSION be;            
			init_BOOLEAN_EXPRESSION(&best, k_num);
			init_BOOLEAN_EXPRESSION(&be, k_num);
			best.mi = 0;

			int t_cnt = 0;
			while (true)        //test all boolean expression
			{
				t_cnt++;
				int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


				for (int k = 0; k < num_sample; k++)
				{
					for (int j = 0; j < k_num; j++)
						be.snp[j] = snp[all_combi[k_num*i+j]].genotype[k] - '0';
					value = be_run(&be, k_num);
					if (1 == value)
					{
						if (k < num_case)
							case_1++;
						else
							control_1++;
					}
					else
					{
						if (k < num_case)
							case_0++;
						else
							control_0++;
					}
					if (case_0 * control_1 > case_1 * control_0)
					{
						be.e_gate = 1;
						be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
					}
					else
					{
						be.e_gate = 0;
						be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
					}
				}
				be.mi = 0;
				be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
				be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
				be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
				be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
				be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
				be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
				if(be.mi > best.mi)
					be_copy(&best, &be, k_num);
				if (1 == next(&be, k_num))
					break;
			}	


			fprintf(ofp, "%d\t%d\t%d\t%d\t%lf\t", i+1, snp[all_combi[k_num*i]].id, snp[all_combi[k_num*i+1]].id,snp[all_combi[k_num*i+2]].id, all_value[i]);
			if (0 < all_value[i])
				fprintf(ofp, "%lf", chisqr(26, 2 * num_sample * LN2 * all_value[i]), combin(num_SNP, k_num));
			else
				fprintf(ofp, "1.0");
			fprintf(ofp, "\t");
			if (0 < all_value[i])
				fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_value[i],26, num_SNP, k_num);
			//fprintf(ofp, "%lf", 2 * num_sample * LN2 * all_value[i]);
			else
				fprintf(ofp, "1.0");

			char s_temp[200];
			make_string(s_temp, best, k_num);
			fprintf(ofp, "%s\n", s_temp); 

			for (int j = 0; j < set_size * k_num; j++)
				for (int k = 0; k < k_num; k++)
					if (temp_set[j] == all_combi[k_num*i+k])
						temp_set[j] = -1;

			be_free(&best);
			be_free(&be);
		}
	}
	fclose(ofp);
	printf("three-locus end\n");    
}

void four_locus()
{
	/*int k_num = 4, opt = 0;

	  printf("############# four-locus #############\n");
	  k_means(k_num, set_size);
	  check_set(k_num);

	  int*temp_set = (int*)malloc(set_size * k_num*sizeof(int));
	  int temp_cnt = 0;

	  for (int i = 0; i < set_size * k_num; i++)
	  {
	  temp_set[i] = centroid[i / set_size].set[i % set_size];
	  if (temp_set[i] >= 0)
	  temp_cnt++;
	  }

	  int* all_combi = (int *)malloc(4*(temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3))
	  / 2 / 3 * sizeof(int));
	  double* all_value = (double *)malloc((temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3))
	  / 2 / 3 /4*sizeof(double));

	  int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3)) / 2 / 3/4;
	  for (int i = 0; i < set_size * k_num - 3; i++)
	  {
	  for (int i1 = i + 1; i1 < set_size * k_num - 2; i1++)
	  for (int i2 = i1 + 1; i2 < set_size * k_num-1; i2++)
	  for (int i3 = i2 + 1; i3 < set_size * k_num; i3++)
	  if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2] && 0 <= temp_set[i3])
	  {
	  all_combi[4*cnt+ 0] = temp_set[i];
	  all_combi[4*cnt+1] = temp_set[i1];
	  all_combi[4*cnt+2] = temp_set[i2];
	  all_combi[4*cnt+3] = temp_set[i3];
	  all_value[cnt] = distance_4(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3]);
	  cnt++;
	  }
	  }
	  printf("four-locus sorting\n");
	//sorting////////////
	for (int i = 0; i < cnt - 1; i++)
	for (int j = i + 1; j < cnt; j++)
	{
	int i_temp;
	double d_temp;
	if (all_value[i] < all_value[j])
	{
	for (int k = 0; k < k_num; k++)
	{ i_temp = all_combi[4*i+k]; all_combi[4*i+k] = all_combi[4*j+k]; all_combi[4*j+k] = i_temp; }
	d_temp = all_value[i]; all_value[i] = all_value[j]; all_value[j] = d_temp;
	}
	}

	FILE *ofp;
	char str_file[200];
	strcpy(str_file, f_case);
	strcat(str_file, "_four_all.txt");

	ofp = fopen(str_file, "w");
	fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tSNP4\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

	for (int i = 0; i < cnt; i++)
	{
	fprintf(ofp, "%d\t%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[4*i]].id, snp[all_combi[4*i+1]].id,
	snp[all_combi[4*i+ 2]].id, snp[all_combi[4*i+ 3]].id, all_value[i]);
	if (0 < all_value[i])
	fprintf(ofp, "%lf", chisqr(80, 2 * num_sample * LN2 * all_value[i])*combin(num_SNP,k_num));
	else
	fprintf(ofp, "1.0");
	fprintf(ofp, "\t");
	if (0 < all_value[i])
	fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_value[i],80, num_SNP, k_num);
	else
		fprintf(ofp, "1.0");

	if (i < best_size)
	{
		BOOLEAN_EXPRESSION best;
		BOOLEAN_EXPRESSION be;
		init_BOOLEAN_EXPRESSION(&best, k_num);
		init_BOOLEAN_EXPRESSION(&be, k_num);
		best.mi = 0;

		int t_cnt = 0;
		while (true)        //test all boolean expression
		{
			t_cnt++;
			int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;

			for (int k = 0; k < num_sample; k++)
			{
				for (int j = 0; j < k_num; j++)
					be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
				value = be_run(&be, k_num);
				if (1 == value)
				{
					if (k < num_case)
						case_1++;
					else
						control_1++;
				}
				else
				{
					if (k < num_case)
						case_0++;
					else
						control_0++;
				}
				if (case_0 * control_1 > case_1 * control_0)
				{
					be.e_gate = 1;
					be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
				}
				else
				{
					be.e_gate = 0;
					be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
				}
			}
			be.mi = 0;
			be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
			be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
			be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
			be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
			be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
			be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
			if (be.mi > best.mi)
				be_copy(&best, &be, k_num);
			if (1 == next(&be, k_num))
				break;

		}
		char s_temp[200];
		make_string(s_temp, best, k_num);
		fprintf(ofp, "%s\n", s_temp); 
		be_free(&best);
		be_free(&be);

	}
	else
		fprintf(ofp,"\n");

}
fclose(ofp);

strcpy(str_file, f_case);
strcat(str_file, "_four_deduplicated.txt");

ofp = fopen(str_file, "w");
fprintf(ofp, "SNP1\tSNP2\tSNP3\tSNP4\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");


for (int i = 0; i < cnt; i++)
{
	int d_cnt = 0;
	for (int j = 0; j < set_size * k_num; j++)
		for (int k = 0; k < k_num; k++)
			if (temp_set[j] == all_combi[k_num*i+ k])
				d_cnt++;

	if (d_cnt == k_num)
	{
		BOOLEAN_EXPRESSION best;
		BOOLEAN_EXPRESSION be;
		init_BOOLEAN_EXPRESSION(&best, k_num);
		init_BOOLEAN_EXPRESSION(&be, k_num);
		best.mi = 0;

		int t_cnt = 0;
		while (true)        //test all boolean expression
		{
			t_cnt++;
			int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


			for (int k = 0; k < num_sample; k++)
			{
				//label12.Text = "for test";
				for (int j = 0; j < k_num; j++)
					be.snp[j] = snp[all_combi[i, j]].genotype[k] - '0';
				value = be_run(&be, k_num);
				//label12.Text = label12.Text + value.ToString();
				//label12.Refresh();
				if (1 == value)
				{
					if (k < num_case)
						case_1++;
					else
						control_1++;
				}
				else
				{
					if (k < num_case)
						case_0++;
					else
						control_0++;
				}
				if (case_0 * control_1 > case_1 * control_0)
				{
					be.e_gate = 1;
					be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
				}
				else
				{
					be.e_gate = 0;
					be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
				}
			}
			be.mi = 0;
			be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
			be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
			be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
			be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
			be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
			be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
			if (be.mi > best.mi)
				be_copy(&best, &be, k_num);
			if (1 == next(&be, k_num))
				break;

		}
		fprintf(ofp, "%d\t%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[4*i]].id, snp[all_combi[4*i+1]].id,
				snp[all_combi[4*i+ 2]].id, snp[all_combi[4*i+ 3]].id, all_value[i]);
		if (0 < all_value[i])
			fprintf(ofp, "%lf", chisqr(80, 2 * num_sample * LN2 * all_value[i])*combin(num_SNP, k_num));
		else
			fprintf(ofp, "1.0");
		fprintf(ofp, "\t");
		if (0 < all_value[i])
			fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_value[i],80, num_SNP, k_num);
		else
			fprintf(ofp, "1.0");

		char s_temp[200];
		make_string(s_temp, best, k_num);
		fprintf(ofp, "%s\n", s_temp); 

		for (int j = 0; j < set_size * k_num; j++)
			for (int k = 0; k < k_num; k++)
				if (temp_set[j] == all_combi[i, k])
					temp_set[j] = -1;

		be_free(&best);
		be_free(&be);
	}
}
fclose(ofp);
printf("four-locus end\n");*/

int k_num = 4;

printf("four-locus\n");
k_means(k_num, best_size);
check_set(k_num);



int* temp_set = (int*)malloc(set_size * k_num*sizeof(int));
int temp_cnt = 0;

for (int i = 0; i < set_size * k_num; i++)
{
	temp_set[i] = centroid[i / set_size].set[i % set_size];
	if (temp_set[i] >= 0)
		temp_cnt++;
}

COMBINATION* all_combi = (COMBINATION*)malloc((temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3)) / 2 / 3 / 4 * sizeof(COMBINATION));

printf("four-locus best_search\n");


int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3)) / 2 / 3 / 4 ;
for (int i = 0; i < set_size * k_num - 3; i++)
{
	for (int i1 = i + 1; i1 < set_size * k_num - 2; i1++)
		for (int i2 = i1 + 1; i2 < set_size * k_num - 1; i2++)
			for (int i3 = i2 + 1; i3 < set_size * k_num; i3++)
				if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2] && 0 <= temp_set[i3])
				{
					all_combi[cnt].snp[0] = temp_set[i];
					all_combi[cnt].snp[1] = temp_set[i1];
					all_combi[cnt].snp[2] = temp_set[i2];
					all_combi[cnt].snp[3] = temp_set[i3];
					all_combi[cnt].value[0] = distance_4(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3]);
					cnt++;
				}
}
printf("four-locus sorting\n");

qsort(all_combi, (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3) ) / 2 / 3 / 4, sizeof(COMBINATION), combiCompare);

FILE *ofp;
char str_file[200];
strcpy(str_file, f_case);
strcat(str_file, "_four_all.txt");

ofp = fopen(str_file, "w");
fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tSNP4\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

for (int i = 0; i < cnt; i++)
{
	fprintf(ofp, "%d\t%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[i].snp[0]].id,
			snp[all_combi[i].snp[1]].id, snp[all_combi[i].snp[2]].id, snp[all_combi[i].snp[3]].id,
			all_combi[i].value[0]);
	if (0 < all_combi[i].value[0])
		fprintf(ofp, "%lf", chisqr(80, 2 * num_sample * LN2 * all_combi[i].value[0])*combin(num_SNP, k_num));
	else
		fprintf(ofp, "1.0");
	fprintf(ofp, "\t");
	if (0 < all_combi[i].value[0])
		fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_combi[i].value[0], 80, num_SNP, k_num);
	else
		fprintf(ofp, "1.0");

	if (i < best_size)
	{
		BOOLEAN_EXPRESSION best;
		BOOLEAN_EXPRESSION be;
		init_BOOLEAN_EXPRESSION(&best, k_num);
		init_BOOLEAN_EXPRESSION(&be, k_num);
		best.mi = 0;

		int t_cnt = 0;
		while (true)        //test all boolean expression
		{
			t_cnt++;
			int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


			for (int k = 0; k < num_sample; k++)
			{
				for (int j = 0; j < k_num; j++)
					be.snp[j] = snp[all_combi[i].snp[j]].genotype[k] - '0';
				value = be_run(&be, k_num);
				if (1 == value)
				{
					if (k < num_case)
						case_1++;
					else
						control_1++;
				}
				else
				{
					if (k < num_case)
						case_0++;
					else
						control_0++;
				}
				if (case_0 * control_1 > case_1 * control_0)
				{
					be.e_gate = 1;
					be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
				}
				else
				{
					be.e_gate = 0;
					be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
				}
			}
			be.mi = 0;
			be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
			be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
			be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
			be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
			be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
			be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
			if (be.mi > best.mi)
				be_copy(&best, &be, k_num);
			if (1 == next(&be, k_num))
				break;

		}
		char s_temp[200];
		make_string(s_temp, best, k_num);
		fprintf(ofp, "%s\n", s_temp);
		be_free(&best);
		be_free(&be);
	}
	else
		fprintf(ofp, "\n");

}
fclose(ofp);

strcpy(str_file, f_case);
strcat(str_file, "_four_deduplicated.txt");

ofp = fopen(str_file, "w");
fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tSNP4\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

for (int i = 0; i < cnt; i++)
{
	int d_cnt = 0;
	for (int j = 0; j < set_size * k_num; j++)
		for (int k = 0; k < k_num; k++)
			if (temp_set[j] == all_combi[i].snp[k])
				d_cnt++;

	if (d_cnt == k_num)
	{
		BOOLEAN_EXPRESSION best;
		BOOLEAN_EXPRESSION be;
		init_BOOLEAN_EXPRESSION(&best, k_num);
		init_BOOLEAN_EXPRESSION(&be, k_num);
		best.mi = 0;

		int t_cnt = 0;
		while (true)        //test all boolean expression
		{
			t_cnt++;
			int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


			for (int k = 0; k < num_sample; k++)
			{
				//label12.Text = "for test";
				for (int j = 0; j < k_num; j++)
					be.snp[j] = snp[all_combi[i].snp[j]].genotype[k] - '0';
				value = be_run(&be, k_num);
				//label12.Text = label12.Text + value.ToString();
				//label12.Refresh();
				if (1 == value)
				{
					if (k < num_case)
						case_1++;
					else
						control_1++;
				}
				else
				{
					if (k < num_case)
						case_0++;
					else
						control_0++;
				}
				if (case_0 * control_1 > case_1 * control_0)
				{
					be.e_gate = 1;
					be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
				}
				else
				{
					be.e_gate = 0;
					be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
				}
			}
			be.mi = 0;
			be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
			be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
			be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
			be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
			be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

			be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
			be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
			if (be.mi > best.mi)
				be_copy(&best, &be, k_num);
			if (1 == next(&be, k_num))
				break;

		}
		fprintf(ofp, "%d\t%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[i].snp[0]].id,
				snp[all_combi[i].snp[1]].id, snp[all_combi[i].snp[2]].id, snp[all_combi[i].snp[3]].id,
				all_combi[i].value[0]);
		if (0 < all_combi[i].value[0])
			fprintf(ofp, "%lf", chisqr(80, 2 * num_sample * LN2 * all_combi[i].value[0])*combin(num_SNP, k_num));
		else
			fprintf(ofp, "1.0");
		fprintf(ofp, "\t");
		if (0 < all_combi[i].value[0])
			fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_combi[i].value[0], 80, num_SNP, k_num);
		else
			fprintf(ofp, "1.0");

		char s_temp[200];
		make_string(s_temp, best, k_num);
		fprintf(ofp, "%s\n", s_temp);

		for (int j = 0; j < set_size * k_num; j++)
			for (int k = 0; k < k_num; k++)
				if (temp_set[j] == all_combi[i].snp[k])
					temp_set[j] = -1;
		be_free(&best);
		be_free(&be);
	}
}
fclose(ofp);
printf("four-locus end\n");
}

void five_locus()
{
	int k_num = 5;

	printf("five-locus\n");
	k_means(k_num, best_size);
	check_set(k_num);



	int* temp_set = (int*)malloc(set_size * k_num*sizeof(int));
	int temp_cnt = 0;

	for (int i = 0; i < set_size * k_num; i++)
	{
		temp_set[i] = centroid[i / set_size].set[i % set_size];
		if (temp_set[i] >= 0)
			temp_cnt++;
	}

	COMBINATION* all_combi = (COMBINATION*)malloc((temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3) * (temp_cnt - 4)) / 2 / 3 / 4 / 5*sizeof(COMBINATION));          

	printf("five-locus best_search\n");


	int cnt = 0, total = (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3) * (temp_cnt - 4)) / 2 / 3 / 4/5;
	for (int i = 0; i < set_size * k_num - 4; i++)
	{
		for (int i1 = i + 1; i1 < set_size * k_num - 3; i1++)
			for (int i2 = i1 + 1; i2 < set_size * k_num - 2; i2++)
				for (int i3 = i2 + 1; i3 < set_size * k_num-1; i3++)
					for (int i4 = i3 + 1; i4 < set_size * k_num; i4++)
						if (0 <= temp_set[i] && 0 <= temp_set[i1] && 0 <= temp_set[i2] && 0 <= temp_set[i3] && 0 <= temp_set[i4])
						{
							all_combi[cnt].snp[0] = temp_set[i];
							all_combi[cnt].snp[1] = temp_set[i1];
							all_combi[cnt].snp[2] = temp_set[i2];
							all_combi[cnt].snp[3] = temp_set[i3];
							all_combi[cnt].snp[4] = temp_set[i4];
							all_combi[cnt].value[0] = distance_5(temp_set[i], temp_set[i1], temp_set[i2], temp_set[i3], temp_set[i4]);                        
							cnt++;
						}
	}
	printf("five-locus sorting\n");

	qsort(all_combi, (temp_cnt * (temp_cnt - 1) * (temp_cnt - 2) * (temp_cnt - 3) * (temp_cnt - 4)) / 2 / 3 / 4/5, sizeof(COMBINATION), combiCompare);    

	FILE *ofp;
	char str_file[200];
	strcpy(str_file, f_case);
	strcat(str_file, "_five_all.txt");

	ofp = fopen(str_file, "w");
	fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tSNP4\tSNP5\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

	for (int i = 0; i < cnt; i++)
	{
		fprintf(ofp, "%d\t%d\t%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[i].snp[0]].id, 
				snp[all_combi[i].snp[1]].id, snp[all_combi[i].snp[2]].id, snp[all_combi[i].snp[3]].id,
				snp[all_combi[i].snp[4]].id, all_combi[i].value[0]);
		if (0 < all_combi[i].value[0])
			fprintf(ofp, "%lf", chisqr(242, 2 * num_sample * LN2 * all_combi[i].value[0])*combin(num_SNP, k_num));
		else
			fprintf(ofp, "1.0");
		fprintf(ofp, "\t");
		if (0 < all_combi[i].value[0])
			fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_combi[i].value[0],242, num_SNP, k_num);
		else
			fprintf(ofp, "1.0");        

		if (i < best_size)
		{
			BOOLEAN_EXPRESSION best;
			BOOLEAN_EXPRESSION be;
			init_BOOLEAN_EXPRESSION(&best, k_num);
			init_BOOLEAN_EXPRESSION(&be, k_num);
			best.mi = 0;

			int t_cnt = 0;
			while (true)        //test all boolean expression
			{
				t_cnt++;
				int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


				for (int k = 0; k < num_sample; k++)
				{
					for (int j = 0; j < k_num; j++)
						be.snp[j] = snp[all_combi[i].snp[j]].genotype[k] - '0';
					value = be_run(&be, k_num);
					if (1 == value)
					{
						if (k < num_case)
							case_1++;
						else
							control_1++;
					}
					else
					{
						if (k < num_case)
							case_0++;
						else
							control_0++;
					}
					if (case_0 * control_1 > case_1 * control_0)
					{
						be.e_gate = 1;
						be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
					}
					else
					{
						be.e_gate = 0;
						be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
					}
				}
				be.mi = 0;
				be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
				be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
				be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
				be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
				be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
				be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
				if (be.mi > best.mi)
					be_copy(&best, &be, k_num);
				if (1 == next(&be, k_num))
					break;

			}                   
			char s_temp[200];
			make_string(s_temp, best, k_num);
			fprintf(ofp, "%s\n", s_temp);   
			be_free(&best);
			be_free(&be);
		}
		else
			fprintf(ofp,"\n");

	}
	fclose(ofp);

	strcpy(str_file, f_case);
	strcat(str_file, "_five_deduplicated.txt");

	ofp = fopen(str_file, "w");
	fprintf(ofp, "rank\tSNP1\tSNP2\tSNP3\tSNP4\tSNP5\tMI\tp-value(in software, after Bonfferoni correction)\tp-value(in excell, after Bonferroni correction)\tBoolean_expression\taccuracy\tsensitivity\tspecificity\tbalanced_accuracy\tTP\tFP\tFN\tTN\n");

	for (int i = 0; i < cnt; i++)
	{
		int d_cnt = 0;
		for (int j = 0; j < set_size * k_num; j++)
			for (int k = 0; k < k_num; k++)
				if (temp_set[j] == all_combi[i].snp[k])
					d_cnt++;

		if (d_cnt == k_num)
		{
			BOOLEAN_EXPRESSION best;
			BOOLEAN_EXPRESSION be;
			init_BOOLEAN_EXPRESSION(&best, k_num);
			init_BOOLEAN_EXPRESSION(&be, k_num);
			best.mi = 0;

			int t_cnt = 0;
			while (true)        //test all boolean expression
			{
				t_cnt++;
				int case_1 = 0, case_0 = 0, control_1 = 0, control_0 = 0, value;


				for (int k = 0; k < num_sample; k++)
				{
					//label12.Text = "for test";
					for (int j = 0; j < k_num; j++)
						be.snp[j] = snp[all_combi[i].snp[j]].genotype[k] - '0';
					value = be_run(&be, k_num);
					//label12.Text = label12.Text + value.ToString();
					//label12.Refresh();
					if (1 == value)
					{
						if (k < num_case)
							case_1++;
						else
							control_1++;
					}
					else
					{
						if (k < num_case)
							case_0++;
						else
							control_0++;
					}
					if (case_0 * control_1 > case_1 * control_0)
					{
						be.e_gate = 1;
						be.TP = case_0; be.TN = control_1; be.FP = control_0; be.FN = case_1;
					}
					else
					{
						be.e_gate = 0;
						be.TP = case_1; be.TN = control_0; be.FP = control_1; be.FN = case_0;
					}
				}
				be.mi = 0;
				be.mi += ((double)case_0 / num_sample) * log(TO_ZERO + (double)case_0 / num_sample) / log(2.0);
				be.mi += ((double)case_1 / num_sample) * log(TO_ZERO + (double)case_1 / num_sample) / log(2.0);
				be.mi += ((double)control_0 / num_sample) * log(TO_ZERO + (double)control_0 / num_sample) / log(2.0);
				be.mi += ((double)control_1 / num_sample) * log(TO_ZERO + (double)control_1 / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + case_1) / num_sample) * log(TO_ZERO + (double)(case_0 + case_1) / num_sample) / log(2.0);
				be.mi += -((double)(control_0 + control_1) / num_sample) * log(TO_ZERO + (double)(control_0 + control_1) / num_sample) / log(2.0);

				be.mi += -((double)(case_0 + control_0) / num_sample) * log(TO_ZERO + (double)(case_0 + control_0) / num_sample) / log(2.0);
				be.mi += -((double)(case_1 + control_1) / num_sample) * log(TO_ZERO + (double)(case_1 + control_1) / num_sample) / log(2.0);
				if (be.mi > best.mi)
					be_copy(&best, &be, k_num);
				if (1 == next(&be, k_num))
					break;

			}
			fprintf(ofp, "%d\t%d\t%d\t%d\t%d\t%d\t%lf\t", (i + 1), snp[all_combi[i].snp[0]].id, 
					snp[all_combi[i].snp[1]].id, snp[all_combi[i].snp[2]].id, snp[all_combi[i].snp[3]].id,
					snp[all_combi[i].snp[4]].id, all_combi[i].value[0]);
			if (0 < all_combi[i].value[0])
				fprintf(ofp, "%lf", chisqr(242, 2 * num_sample * LN2 * all_combi[i].value[0])*combin(num_SNP, k_num));
			else
				fprintf(ofp, "1.0");
			fprintf(ofp, "\t");
			if (0 < all_combi[i].value[0])
				fprintf(ofp, "=chisq.dist.rt(%lf,%d)*combin(%d,%d)", 2 * num_sample * LN2 * all_combi[i].value[0],242, num_SNP, k_num);
			else
				fprintf(ofp, "1.0");    

			char s_temp[200];
			make_string(s_temp, best, k_num);
			fprintf(ofp, "%s\n", s_temp); 

			for (int j = 0; j < set_size * k_num; j++)
				for (int k = 0; k < k_num; k++)
					if (temp_set[j] == all_combi[i].snp[k])
						temp_set[j] = -1;
			be_free(&best);
			be_free(&be);
		}
	}
	fclose(ofp);
	printf("five-locus end\n");
}

int main(void){
	snp = (SNP*)malloc(SIZE*sizeof(SNP));

	loadData();
	if(from<= 1 && 1 <= to) 
		single_locus();
	if(from<= 2 && 2 <= to) 
		two_locus();
	if(from<= 3 && 3 <= to) 
		three_locus();
	if(from<= 4 && 4 <= to) 
		four_locus();
	if(from<= 5 && 5 <= to) 
		five_locus();

	return 0;
}

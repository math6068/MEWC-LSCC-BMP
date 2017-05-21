// =====================================================================================
//
//       Filename:  LSCC+BMP.cpp
//
//    Description:  This is a solver for weighted maximum clique problem based on SCC and BMP
//
//        Version:  1.0
//        Created:  2016.01
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Yiyuan Wang          yiyuanwangjlu@126.com
//                  Shaowei Cai          caisw@ios.ac.cn
//                  Minghao Yin          ymh@nenu.edu.cn

//         For example：./LSCC+BMP bio-yeast.mtx 1000 1 

// *************************************************************
//         Date: 2016.10
//     Revision: Extend LSCC+BMP from handling maximum vertec weight clique problem to handling maximum edge weight clique problem.
//       Author: Yi Fan   yi.fan4@griffithuni.edu.au     
// =====================================================================================

      
      
#include<limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include<unistd.h>
#include <time.h>
#include <ctime>
//#include <vector>
#include <string.h>
#include<math.h>
#include <climits>
#include <float.h>
#include<assert.h>

#define linux_mode

#ifdef linux_mode
#include <sys/times.h>
#endif
#ifdef windows_mode
//#include <windows.h>
#endif

using namespace std;

typedef long long LL;

//const long double EPS = 0.0001; //10^-4

const long double EPS = 0.000000000001; //10^-12


// swith between BMS and noBMS mode
//int BMP = INT_MAX;
int BMP = 100;


#ifdef linux_mode
tms start, finish;
double time_limit;
double real;
#endif
#ifdef windows_mode
int time_limit;
int real;
unsigned long start_time;
#endif


//int lbest;
long double lbest;
     int M_iter = 0;
struct Edge1{
	int v1;
	int v2;
};

Edge1 *edge;  //MAXE
int *v_degree_tmp;//MAXV

int *adjaclen;//MAXV

//double real_solve1=-1;
unsigned long real_solve1=-1;
//double real_solve2=-1;
unsigned long real_solve2=-1;
int **neighbor;
int *neighbor_len;//MAXV
int **adj_edge;
int *conf_change;//MAXV

int *time_stamp;//MAXV




int *temp_array;//MAXV

char * File_Name;

int *vectex;//MAXV
int *funch;//MAXV
int *address;//MAXV
int *tabuin;//MAXV
int Max_Vtx,  Max_Iter; 
//int f;
//int fbest;

int *cruset;//MAXV
int len;
int tm1;
int tm2;
int *C0; //MAXV// the candidates node for ADD operation?
int *C1; //MAXV// the Candidate nodes for SWAP operation?
//int *We; //MAXV// Weights of nodes or edges
long double *We;
long double *score;
int *BC;//MAXV // 
int len0; // the length of C0
int len1; // the length of C1
int *TC1;//MAXV // Temporal candidate nodes?
int Iter; // the number of iterations taken
int TABUL = 7;
//int Wf;
long double Wf;
//int Wbest;
long double Wbest;
int *FC1;//MAXV
int *Tbest;//MAXV
int *TTbest;//MAXV
//int Waim;
long double Waim;
int Titer;

int len_best = 0;
int len_W;
int iter_best;
int Iteration[ 100 ];
double time_used[ 100 ];
int len_used[ 100 ];
//int W_used[ 100 ];
char outfilename[30];
int len_improve;
int len_time;
//int Wmode;
//int TABUL0 = 5;

/************************************************************************/
/*   WYY: Variables for Configuration Checking                    */
/************************************************************************/
/* neighbor[i][j] = n means that i and n are connceted, and n is the jth 
 *  neighbor of n.
 */

/* time_stamp[m]=j means the conf_change value of node m recently 
 * changes at jth iteration.
 */

//struct  timeval start;
//struct  timeval end;

void dump_neighborhood();
void show_state();
int verify();
int validate();
#define DEBUG 0
///////////////////////////

inline int long_double_equals(long double a, long double b, long double epsilon = EPS)
{
	return fabs(a - b) < epsilon;
}

inline int long_double_greater(long double a, long double b, long double epsilon = EPS)
{
	return a - b >= epsilon;
}

inline int long_double_smaller(long double a, long double b, long double epsilon = EPS)
{
	return b - a >= epsilon;
}

inline int long_double_geq(long double a, long double b, long double epsilon = EPS)
{
	return  (a - b >= epsilon || fabs(a - b) < epsilon);
}

inline int long_double_leq(long double a, long double b, long double epsilon = EPS)
{
	return (b - a >= epsilon || fabs(a - b) < epsilon);
}

inline int edge_is(int m, int n)
{
  int i;
  int index;
  for(i=0;i<neighbor_len[m];i++){
    index=neighbor[m][i];
    if(index==n)return 0;
  }
  return 1;
}


//int **adj_mat, int *v_weight, int _v_num, int _e_num
// section 0, initiaze

void Initializing()
{
     ifstream FIC;
     FILE *fp;
     FIC.open(File_Name);
     if ( FIC.fail() )
     {
           cout << "### Erreur open, File_Name " << File_Name << endl;
           getchar();
           exit(0);
     }
     char line[1024];
     //Max_Vclique=300;

     if ( FIC.eof() )
     {
           cout << "### Error open, File_Name " << File_Name << endl;
           exit(0);
     }
     int nb_vtx=0, nb_edg=-1, max_edg=0;
     int x1, x2;
	 int tt; 

    //FIC >> StrReading;

/*
	do{
	FIC.getline(StrReading, 100);
	}while(StrReading[0] != 'p');
	char tmpstr1[10], tmpstr2[10];
	sscanf(StrReading, "%s %s %d %d", tmpstr1, tmpstr2, &Max_Vtx, &nb_edg);
*/
	FIC.getline(line, 1024);
	sscanf(line, "%d %d", &Max_Vtx, &nb_edg);
//cout << Max_Vtx << " " << nb_edg << endl;
	//}while(tt != Max_Vtx);
	 //FIC >> tt >> Max_Vtx >> nb_edg;
	 //neighbor=(int **)malloc(nb_edg*2*sizeof(int));//neighbor set
	//neighbor = new int*[nb_edg * 2];
//cout << 1 << endl;
		neighbor = new int*[Max_Vtx];
		adj_edge = new int*[Max_Vtx];
//cout << 2 << endl;
/////////////////////////////////////
//cout << nb_edg << endl;
	edge = new Edge1[nb_edg];
	We = new long double[nb_edg];
	score = new long double[Max_Vtx];
//cout << 2.5 << endl;
	v_degree_tmp = new int[Max_Vtx];
	adjaclen = new int[Max_Vtx];//MAXV
//cout << 3 << endl;
	neighbor_len = new int[Max_Vtx];//MAXV
	conf_change = new int[Max_Vtx];//MAXV
	time_stamp = new int[Max_Vtx];//MAXV
//cout << 4 << endl;
	temp_array = new int[Max_Vtx];//MAXV
	vectex = new int[Max_Vtx];//MAXV
	funch = new int[Max_Vtx];//MAXV
//cout << 5 << endl;
	address = new int[Max_Vtx];//MAXV	
	tabuin = new int[Max_Vtx];//MAXV	
	cruset = new int[Max_Vtx];//MAXV	
	C0 = new int[Max_Vtx];//MAXV	
	C1 = new int[Max_Vtx];//MAXV	
	//We = new int[Max_Vtx];//MAXV
	//We = new long double[Max_Vtx];//MAXV	
	BC = new int[Max_Vtx];//MAXV	
	TC1 = new int[Max_Vtx];//MAXV	
	FC1 = new int[Max_Vtx];//MAXV	
	Tbest = new int[Max_Vtx];//MAXV	
	TTbest = new int[Max_Vtx];//MAXV	
/////////////////////////////////////

	 for( int x = 0 ; x < Max_Vtx ; x++ ) 
	 {
		 conf_change[x] = 1; // initialize
         time_stamp[x] = 0;
		 adjaclen[x]=0;
         neighbor_len[x]=0;
         vectex[x]  = x;
         address[x] = x;
	}
	//char tmpstr[10];
	int tmpint;
	int ii;
/*
	for(ii = 0; ii < Max_Vtx; ii++)
	{
		int v_index;
		//int v_weight;
		long double v_weight;
		FIC >> tmpstr >> v_index >> v_weight;
		We[v_index - 1] = v_weight;
//cout << v_index << " " << v_weight << endl;
//getchar();
	}
*/
	for(ii=0;ii<nb_edg;ii++)
	{
		long double e_weight;

		FIC >> tmpint >> x1 >> x2 >> e_weight;

		x1--; x2--;
        if ( x1<0 || x2<0 || x1>=Max_Vtx || x2 >=Max_Vtx )
        {
            cout << "### Error of node : x1="
                 << x1 << ", x2=" << x2 << endl;
            exit(0);
        }
				if(tmpint != ii)
				{
						cout << "### Error edge: " << x1 << " and " << x2 << endl;
						exit(0);
				}
		max_edg++;
		neighbor_len[x1]++;
		neighbor_len[x2]++;
		edge[ii].v1 = x2;
		edge[ii].v2 = x1;
		We[ii] = e_weight;     
		}
    int v;
		for	(v=0; v<Max_Vtx; v++) score[v]= 0;
    for (v=0; v<Max_Vtx; v++) adjaclen[v]=Max_Vtx-1-neighbor_len[v];//????

	for (v=0; v<Max_Vtx; v++)
	{ 
		//neighbor[v]=(int *)malloc( neighbor_len[v]*sizeof(int));
		neighbor[v]= new int[neighbor_len[v]];
		adj_edge[v]= new int[neighbor_len[v]];
	}

    for(v=0; v<Max_Vtx; v++) v_degree_tmp[v]=0; 
	int e,v1,v2;  
	for (e=0; e<nb_edg; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;
		neighbor[v1][v_degree_tmp[v1]] = v2;
		neighbor[v2][v_degree_tmp[v2]] = v1;
		adj_edge[v1][v_degree_tmp[v1]] = e;
		adj_edge[v2][v_degree_tmp[v2]] = e;
		v_degree_tmp[v1]++;
		v_degree_tmp[v2]++;
	}
	int i;


     if ( 0 && max_edg != nb_edg )
     {
           cout << "### Error de lecture du graphe, nbre aretes : annonce="
                 << nb_edg << ", lu=" << max_edg  << endl;
           exit(0);
     }
     
     

     for( int x = 0; x < Max_Vtx; x++ )
     {
        //We[ x ] = (x+1)%Wmode + 1;
        BC[ x ] = 0;
        //We[ x ] = 1;
        //We[ x ] = ( rand() % 10 ) + 1;
     }
    
     FIC.close();     
}

/*
void Initializing(int **adj_mat, long double *v_weight, int _v_num, int _e_num)
{
     int nb_vtx=0, nb_edg=-1, max_edg=0;
     int x1, x2;
	 int tt; 

	Max_Vtx = _v_num;
	nb_edg = _e_num;

		neighbor = new int*[Max_Vtx];


	edge = new Edge1[nb_edg];

	v_degree_tmp = new int[Max_Vtx];
	adjaclen = new int[Max_Vtx];//MAXV

	neighbor_len = new int[Max_Vtx];//MAXV

	conf_change = new int[Max_Vtx];//MAXV
	time_stamp = new int[Max_Vtx];//MAXV

	temp_array = new int[Max_Vtx];//MAXV

	vectex = new int[Max_Vtx];//MAXV
	funch = new int[Max_Vtx];//MAXV

	address = new int[Max_Vtx];//MAXV	

	tabuin = new int[Max_Vtx];//MAXV	
	cruset = new int[Max_Vtx];//MAXV	

	C0 = new int[Max_Vtx];//MAXV	
	C1 = new int[Max_Vtx];//MAXV	

	//We = new int[Max_Vtx];//MAXV
	We = new long double[Max_Vtx];//MAXV	

	BC = new int[Max_Vtx];//MAXV	
	TC1 = new int[Max_Vtx];//MAXV	

	FC1 = new int[Max_Vtx];//MAXV	
	Tbest = new int[Max_Vtx];//MAXV	

	TTbest = new int[Max_Vtx];//MAXV	
/////////////////////////////////////


	 for( int x = 0 ; x < Max_Vtx ; x++ ) 
	 {
		 conf_change[x] = 1; // initialize

         time_stamp[x] = 0;
		 adjaclen[x]=0;

         neighbor_len[x]=0;
         vectex[x]  = x;

         address[x] = x;
	}

	char tmpstr[10];

	for(int i = 1; i <= v_num; i++)
	{
		We[i - 1] = v_weight[i];
	}


	int e = 0;
	for(int i = 1; i <= v_num; i++)
		for(int j = i + 1; j <= v_num; j++)
		{
			if(adj_mat[i][j] == 0) continue;

			v1 = i - 1; v2 = j - 1;

			neighbor_len[v1]++;
			neighbor_len[v2]++;

			edge[e].v1 = v1;
			edge[e].v2 = v2;  

			e++;
		}	

    int v;

    for (v=0; v<Max_Vtx; v++) adjaclen[v]=Max_Vtx-1-neighbor_len[v];//????


	for (v=0; v<Max_Vtx; v++) 
		neighbor[v]= new int[neighbor_len[v]];

    for(v=0; v<Max_Vtx; v++) v_degree_tmp[v] = 0; 
	int v1,v2;  
	for (e=0; e<nb_edg; e++)
	{

		v1=edge[e].v1;
		v2=edge[e].v2;

		neighbor[v1][v_degree_tmp[v1]] = v2;
		neighbor[v2][v_degree_tmp[v2]] = v1;

		v_degree_tmp[v1]++;
		v_degree_tmp[v2]++;

	}
	int i;

     if ( 0 && max_edg != nb_edg )
     {
           cout << "### Error de lecture du graphe, nbre aretes : annonce="
                 << nb_edg << ", lu=" << max_edg  << endl;

           exit(0);
     }

     for( int x = 0; x < Max_Vtx; x++ )
     {

        BC[ x ] = 0;

     }

}
*/
void free_memory()
{
	for(int v=0; v<Max_Vtx; v++)
	{
		delete[] neighbor[v];
		delete[] adj_edge[v];
	}
	delete[] neighbor;
	delete[] adj_edge;
	delete[] edge;
	delete[] v_degree_tmp;
	delete[] adjaclen;
	delete[] neighbor_len;
	delete[] conf_change;
	delete[] time_stamp;
	delete[]  temp_array;
	delete[] vectex;
	delete[]  funch;
	delete[] address;
	delete[] tabuin;
	delete[] cruset;
	delete[] C0;
	delete[] C1;
	delete[] We;
	delete[] BC;
	delete[] TC1;
	delete[] FC1;
	delete[] Tbest;
	delete[]  TTbest;
	delete[] score;
}

// WYY
void dump_conf_change() {
	printf("\nconf_change:\n");
	for(int i = 0; i < Max_Vtx; i++) {
		printf("%4d(%d) ", i, conf_change[i]);
	}
	printf("\n");
}

// WYY
inline void neighbor_add(int node) {
	int node2;
	int num_neighbor = neighbor_len[node];
	
	conf_change[node] = 0; // keep this node from being removed immediately
	time_stamp[node] = Iter;
	for(int i = 0; i < num_neighbor; i++) {
		node2 = neighbor[node][i];
		if(conf_change[node2] == 0){
			conf_change[node2] = 1;
			time_stamp[node2] = Iter;//????
		}
	}
}
// WYY
inline void neighbor_drop(int node) {
	conf_change[node] = 0;
	time_stamp[node] = Iter;
}

// WYY
inline bool is_forbiden_cc(int node) {
	return (conf_change[node] == 0 ? true : false);
}


// WYY
void dump_neighborhood() {
	printf("Neighborhood:\n");
	for(int i = 0; i < Max_Vtx; i++) {
		printf(": ");
		for(int j = 0; j < neighbor_len[i]; j++)
			printf("%d ", neighbor[i][j]);
		printf("\n");	
	}
	return;
}

// WYY
void dump_cur_clique() {
	return;
	int n;
	printf("\ncurrent clique includes %d nodes:", len);
	for(int i = 0; i < len; i++) {
		n = cruset[i];
		printf("%d(%d) ", n, vectex[n]);
	}
	printf("\n");
}

int randomInt( int n )
{
    return rand() % n;
}

void clearGamma()
{
    int i, j, k, l;
    memset( vectex, 0, tm1 );
    memset(  funch, 0, tm1 );
    memset(address, 0, tm1 );
    memset( tabuin, 0, tm1 );
		memset(score, 0, sizeof(long double) * Max_Vtx);
    for( i = 0; i < Max_Vtx; i++ )
    {
       C0[ i ] = i;// at the beginning all vertices can be added
       address[ i ] = i;
    }
    len0 = Max_Vtx;
    len1 = 0;
    len = 0;
    Wf = 0;
    Wbest = 0;
}

// WYY: C0 is the set of nodes that can be added? and C1 is the set nodes that can be swapped with?
int selectC0( ) // select a random vertex from add_set s.t. confChange=1, if called right after intialization, the confChange info is ignored
//breaking ties randomly
{
    int i, j, k, l, m;
    l = 0;
    
    if( len0 > 30 )
    {
       k = randomInt( len0 );
       return k;
    }
    
    // WYY: TC1 records the set of nodes which are not being forbidden
    for( i = 0; i < len0; i++ )
    {
       k = C0[ i ];
       if( !is_forbiden_cc(k) ) // Added by: WYY
       {
         TC1[ l++ ] = i;
       }
    }
    
    if( l == 0 )
      return -1;
    else
    {
        k = randomInt( l );
        k = TC1[ k ];
        return k;
    }
}
// WYY: Select from C0, a node, which is not in tabu and has the max weight, 
// or satisfy the aspiration rule though in tabu.
int WselectC0( )//breaking ties in favor of the greatest age
{
    int i, j, k, l1, l2, m;
		long double score1, score2;
    l1 = 0;
    l2 = 0;
    score1 = 0;
    score2 = 0;
    
    for( i = 0; i < len0; i++ )
    {
       k = C0[ i ];
       
       // WYY:store nodes that are not in tabu list and with the maximum weight, in FC1
       if( !is_forbiden_cc(k) ) // Added by WYY
       {
           //if( We[ k ] > w1 )
			if(long_double_greater(score[k], score1))
           {
              l1 = 0;
              score1 = score[ k ];
              FC1[ l1++ ] = i;
           }
           //else if ( We[ k ] >= w1 )
			else if(long_double_geq(score[k], score1))
           {
              FC1[ l1++ ] = i;// free best
           }
       }
       else
       {  // WYY: stores nodes that are being in tabu but with the maximum weight, in TC1
          //if( We[ k ] > w2 )
			if(long_double_greater(score[k], score2))
           {
              l2 = 0;
              score2 = score[ k ];
              TC1[ l2++ ] = i;
           }
           //else if ( We[ k ] >= w2 )
			else if(long_double_geq(score[k], score2))
           {
              TC1[ l2++ ] = i;// tabu best
           }
       }
    }
    
    // WYY: to check first if the aspiration rule is applicable.
    // If not, select a nodes which have the highest weithgts; break ties randomly.
    //if( (l2 > 0) && ( w2 > w1 ) && ((w2+Wf)>Wbest) ) //tabu best can lead to a new best clique
	if( l2 > 0 && long_double_greater(score2, score1) && long_double_greater((score2+Wf), Wbest) )
    {
        /*
        k = randomInt( l2 );
        k = TC1[ k ];
        */
        
        // WYY: Select the node with the oldest age
        k = TC1[0];
        int oldest_time = time_stamp[ C0[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l2; j++) {
			index = TC1[j];
			node = C0[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
			}
		}
		      
        //cout << "yes in aspiration w2+Wf = " << w2+Wf << endl;
        //getchar();
        return k;
    }  
    //else if( l1 > 0 ) //select free best
	else if(l1 > 0)
    {
    	/*
        k = randomInt( l1 );
        k = FC1[ k ];
        */
        
        // WYY
        k = FC1[0];
        int oldest_time = time_stamp[ C0[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l1; j++) {
			index = FC1[j];
			node = C0[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
				//cout << "elder one found" << endl;
				//getchar();
			}
		}
        return k;
    }
    else
    {
        return -1;
    }
}

// SelN: the index of the node selected from C0
int expand(int SelN)
{

    int i, j, k, k1, l, am, m, n, n1;
    
    m = C0[ SelN ]; // the node is m
//cout << "expand with: " << m + 1 << endl;
    cruset[ len++ ] = m; // add m into the current set
    vectex[ m ] = 1; // set the flag?
    Wf = Wf + score[ m ]; // Wf is the weight of the current clique, i,e, weight found. update it.

		score[m] = -score[m];

	 /* WYY: set the nodes in cruset as neighbors of m, and vice versa */
	if(DEBUG) {
		printf("\nin expand");
		printf("\nadd node %d", m);
	}
	//neighbor_add(m);
	// remove from add-set    
    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;
    
    
    
    int node2;	
    for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    
    	conf_change[m] = 0; // keep this node from being removed immediately
	time_stamp[m] = Iter;
    for(i=0;i<neighbor_len[m];i++){
     node2 = neighbor[m][i];
			int adj_e = adj_edge[m][i];
			if(vectex[node2])
			{
				score[node2] -= We[adj_e];
//cout << "score[" << node2 + 1 << "] -= " << We[adj_e] << endl;
			}
			else
			{
				score[node2] += We[adj_e];
//cout << "score[" << node2 + 1 << "] += " << We[adj_e] << endl;
			}
     if(conf_change[node2] == 0){
		conf_change[node2] = 1;
		time_stamp[node2] = Iter;// time stamp records at what time configration changes, different from the one described in the paper
	}
      n=neighbor[m][i];
      temp_array[n]=1;// temp_array[n]=1 means neighbors of the added vertex
    }

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m)continue;
       if(temp_array[i]==1)continue;
       n=i;
       funch[ n ]++; // WYY: funch[n] traces the number of nodes that are in the current clique
       
       if( funch[ n ] == 1 )
       {   // WYY: remove n from C0 (add-set)
           k1 = address[ n ];
           len0--;
           n1 = C0[ len0 ];
           C0[ k1 ] = n1;
           address[ n1 ] = k1;
           
           // put it into C1 (swap-set)
           C1[ len1 ] = n;
           address[ n ] = len1;
           len1++;
           
           BC[ n ] = m; // WYY: BC[n] = m denotes that n is Being Connected by m.// n is being paired with m
       }
       else if( funch[ n ] == 2 )
       {
           // remove n it from C1
           len1--;
           n1 = C1[ len1 ];
           k1 = address[ n ];
           C1[ k1 ] = n1;
           address[ n1 ] = k1;
       }
    }
/*
cout << "Wf: " << Wf << endl;
show_state();
if(!validate()) exit(1);    
*/
    //if( Wf > Wbest )
	if(long_double_greater(Wf, Wbest))
     {
#ifdef linux_mode
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
				real_solve2 = timeGetTime() - start_time;
#endif
        real_solve2 = round(real_solve2 * 100)/100.0; 
//cout << "Wf: " << Wf << endl;
        Wbest = Wf;
        len_best = len;
//cout << "Wbest: " << Wbest << endl;
		//if(Wbest > lbest)
		if(long_double_greater(Wbest, lbest))
		{
			iter_best = M_iter + Iter;
//cout << "update best solution" << endl;
        for( i = 0; i < Max_Vtx; i++ )
        {
            Tbest[ i ] = vectex[ i ];
        }

     }
		}
    
    return 1;   
}

int selectC1( )// randomly select a confChanged vertex
{
    int i, j, k, l, m;
    l = 0;
    for( i = 0; i < len1; i++ )
    {
       k = C1[ i ];
       if( !is_forbiden_cc(k) ) // WYY
       {
         TC1[ l++ ] = i;
       }
    }
    if( l == 0 )
      return -1;
    else
    {
        k = randomInt( l );
        k = TC1[ k ];
        return k;
    }
}

int kkk;
int WselectC1( )
{
     int i, j, k, l, l1, l2, m, n;
			long double wmn, score1, score2;
     l1 = 0;
     l2 = 0;
     score1 = -1000000;
     score2 = -1000000;
     l = 0;
     for( i = 0; i < len1; i++ )// for each pair in the swap-set
     {
         m = C1[ i ];
         n = BC[ m ];
         if( (vectex[ n ] == 1) && (edge_is( m, n)==1) )// vectex[n]==1 means in clique
           l++;
         else
         {
             for( j = 0; j < len; j++ )
             {
                k = cruset[ j ];
                if( edge_is(m, k)== 1 )
                  break;
             }
             BC[ m ] = k;
         }
     }//maintain BC
  //wyy-20150601
     int count;
     kkk=BMP;
	 if(len1<=BMP) kkk=len1;
	 
	 if(len1==0)return -1;
     //---add end
     for( count = 0; count < kkk; count++ )
     {
		 if(len1>BMP) i=rand()%len1;
		 else i=count;
		 
         m = C1[ i ];
         n = BC[ m ];
         //wmn = We[ m ] - We[ n ];
				wmn = score[m] + score[n];
		 if( !is_forbiden_cc(m) ) // WYY//tabu best
         { // find the nodes that lead highest weight-increase for the current clique.
             //if( wmn > w1 )
			if(long_double_greater(wmn, score1))
             {
                l1 = 0;
                score1 = wmn;
                FC1[ l1++ ] = i;
             }
             //else if ( wmn >= w1 )
			else if(long_double_geq(wmn, score1))
             {
                FC1[ l1++ ] = i;
             }
         }
         else//free best
         {
             //if( wmn > w2 )
			if(long_double_greater(wmn, score2))
             {
                l2 = 0;
                score2 = wmn;
                TC1[ l2++ ] = i;
             }
             //else if ( wmn >= w2 )
			else if(long_double_geq(wmn, score2))
             {
                TC1[ l2++ ] = i;
             }
         }
     }
     
     //if( (l2 > 0) && ( w2 > w1 ) && ((w2+Wf)>Wbest) )//free best
	if(l2 > 0 && long_double_greater( score2, score1 ) && long_double_greater((score2+Wf), Wbest))
     {
        /* 
        k = randomInt( l2 );
        k = TC1[ k ];
        */
        
        // WYY: Select the oldest node
        k = TC1[0];
        int oldest_time = time_stamp[ C1[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l2; j++) {
			index = TC1[j];
			node = C1[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
			}
		}
        return k;
     }  
     //else if( l1 > 0 )//tabu best
	else if(l1 > 0)
     {
     /*
        k = randomInt( l1 );
        k = FC1[ k ];
       */  
        
        // WYY: Select the oldest node
        k = FC1[0];
        int oldest_time = time_stamp[ C1[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l1; j++) {
			index = FC1[j];
			node = C1[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
				//cout << "elder nodes found " << endl;
				//getchar();
			}
		}
       	
       	return k;
     }
     else
     {
         return -1;
     }
}

int plateau( int SelN )
{
     int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;
     
     m = C1[ SelN  ];
     // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
     for(ti = 0; ti < len; ti++)
     {
         m1 = cruset[ ti ];
         if( edge_is(m1, m)== 1 )
            break;
     }
     
     Wf = Wf + score[ m ] + score[ m1 ];
			score[m] = -score[m];
			score[m1] = -score[m1];
/*
cout << "plateau with: " << m + 1 << " and " << m1 + 1 << endl;
cout << "Wf: " << Wf << endl;
show_state();
*/
     //the expand process, put m into the current independent set
     vectex[ m ] = 1;
     cruset[ len++ ] = m;

	 /* WYY: set the nodes in cruset as neighbors of m, and vice versa */
	 if(DEBUG) {
	 	printf("\nin plateau: add node %d", m);
	 	dump_cur_clique();
	 }
	 //neighbor_add(m); // Attention: here, we don't change conf_change values for m's neighbors
	 
     //delete m from C1
     k1 = address[ m ];
     len1--;
     n1 = C1[ len1 ];
     C1[ k1 ] = n1;
     address[ n1 ] = k1;
     
     
    for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    for(i=0;i<neighbor_len[m];i++){
     
      n=neighbor[m][i];
			int e=adj_edge[m][i];
      temp_array[n]=1;
			if(vectex[n])
			{
				score[n] -= We[e];
			}
			else
			{
				score[n] += We[e];
			}
    }
/*
cout << "after adding: " << m + 1 << endl;
show_state();
*/
    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m)continue;
       if(temp_array[i]==1)continue;
       n=i;
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )//vectex[n]=0 means that n is not in clique
        {
             //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
             k1 = address[ n ];
             len0--;
             n1 = C0[ len0 ];
             C0[ k1 ] = n1;
             address[ n1 ] = k1;//remove from add-set
             
             C1[ len1 ] = n;
             address[ n ] = len1;
             len1++;
             BC[ n ] = m;//add into swap-set
           
             //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;//remove from swap-set
        }        
     } 
     
     //the backtrack process, delete m1 from the current independent set
     vectex[ m1 ] = 0;

	 len--;
     cruset[ ti ] = cruset[ len ];
     C1[ len1 ] = m1;
     address[ m1 ] = len1;
     len1++;
     
     /* WYY: neighborhood updating */
	 if(DEBUG) {
	 	printf("\nin plateau: remove node %d", m1);
	 	dump_neighborhood();
	 }
	 neighbor_drop(m1);
	 
    for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    for(i=0;i<neighbor_len[m1];i++){
     
      n=neighbor[m1][i];
			int e=adj_edge[m1][i];
      temp_array[n]=1;
			if(vectex[n])
			{
				score[n] += We[e];
			}
			else
			{
				score[n] -= We[e];
			}
    }//compute the neighbor set to exclude

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m1)continue;
       if(temp_array[i]==1)continue;
       n=i;
	 
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
           k1 = address[ n ];           
           len1--;
           n1 = C1[ len1 ];
           C1[ k1 ] = n1;
           address[ n1 ] = k1;//remove from swap-set
           
           C0[ len0 ] = n;
           address[ n ] = len0;
           len0++;//add into add-set
        }
        else if( funch[ n ] == 1 )//add into swap-set
        {
           C1[ len1 ] = n;
           address[ n ] = len1;
           len1++;
        }
     }
/*
cout << "after dropping: " << m1 + 1 << endl;
cout << "Wf: " << Wf << endl;
show_state();
if(!validate()) exit(1);
*/
/*
show_state(); 
*/  
     //if( Wf > Wbest )
	if(long_double_greater(Wf, Wbest))
     {
#ifdef linux_mode
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
				real_solve2 = timeGetTime() - start_time;
#endif
        real_solve2 = round(real_solve2 * 100)/100.0; 
//cout << "Wf: " << Wf << endl;
        Wbest = Wf;
        len_best = len;
//cout << "Wbest: " << Wbest << endl;
		//if(Wbest > lbest)
		if(long_double_greater(Wbest, lbest))
		{
			iter_best = M_iter + Iter;
//cout << "update best solution" << endl;
        for( i = 0; i < Max_Vtx; i++ )
        {
            Tbest[ i ] = vectex[ i ];
        }
		}
     }
     
     
     return 1;   
}

// WYY: find nodes with minimum weight in the current clique to remove
int Mumi_Weigt()
{
    int i, j, k, l1, m;
    long double w1 = 5000000;
    l1 = 0;
    // WYY: find in cruset the nodes with lowest weights, breaking ties in favor of the oldest one
    for( i = 0; i < len; i++ )
    {
       k = cruset[ i ];
       //if( We[ k ] < w1 )
		if(long_double_smaller(We[k], w1))
       {
          l1 = 0;
          w1 = We[ k ];
          FC1[ l1++ ] = i;
       }
       //else if ( We[ k ] <= w1 )
	else if ( long_double_leq(We[ k ], w1))
       {
          FC1[ l1++ ] = i;
       }
    }
    
    if( l1 == 0 )
      return -1;
    /*
    k = randomInt( l1 );
    k = FC1[ k ];
    */
    
    // WYY: remove an oldest node
    k = FC1[0];
    int oldest_time = time_stamp[ cruset[k] ];
    int cur_index;
    int cur_node;
    int cur_time;
    for(int i = 1; i < l1; i++) {
    	cur_index = FC1[i];
    	cur_node = cruset[cur_index ];
    	cur_time = time_stamp[cur_node];
    	if(cur_time < oldest_time) {
    		oldest_time = cur_time;
    		k = cur_index;
    		//cout << "elder node found " << endl;
    		//getchar();
    	}
    }
    return k;
}

int backtract()
{
     int i, j, k, l, m, m1, n, ti, k1, n1;
     ti = Mumi_Weigt();
     if( ti == -1 )
      return -1;
      
     if(DEBUG) printf("in backtrack");
     
     m1 = cruset[ ti ];

     Wf = Wf + score[ m1 ];

/*
cout << "score: " << score[m1] << endl;
cout << "2. Wf: " << Wf << endl;
getchar();
*/
		 score[m1] = -score[m1];
     vectex[ m1 ] = 0;

     len--;
     cruset[ ti ] = cruset[ len ];
/*
cout << "backtract with: " << m1 + 1 << endl;
cout << "Wf: " << Wf << endl;
if(!validate()) exit(1);
*/     
     /* WYY: functions of neighborhood updating */
	 neighbor_drop(m1);
	 if(DEBUG) {
	 	printf("\nremove node %d", m1);
	 	dump_cur_clique();
	 }
	 //dump_conf_change();
	 //getchar();
     
     C0[ len0 ] = m1;
     address[ m1 ] = len0;
     len0++;
     
   for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    for(i=0;i<neighbor_len[m1];i++){
     
      n=neighbor[m1][i];
			int e=adj_edge[m1][i];
      temp_array[n]=1;
			if(vectex[n])
			{
				score[n] += We[e];
			}
			else
			{
				score[n] -= We[e];
			}
    }

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m1)continue;
       if(temp_array[i]==1)continue;
       n=i;     

        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
           k1 = address[ n ];           
           len1--;
           n1 = C1[ len1 ];
           C1[ k1 ] = n1;
           address[ n1 ] = k1;//remove from swap-set
           
           C0[ len0 ] = n;
           address[ n ] = len0;
           len0++;//add into add-set
        }
        else if( funch[ n ] == 1 )
        {
           C1[ len1 ] = n;
           address[ n ] = len1;
           len1++;// add into add-set
        }
     }
}

void show_state()
{
	cout << "current clique: " << endl;
	for(int v = 0; v < Max_Vtx; v++)
	{
		if(vectex[v]) cout << v + 1 << "\t";
	}
	cout << endl;
	cout << "clique weight: " << Wf << endl;
	cout << "score: " << endl;
	for(int v = 0; v < Max_Vtx; v++)
	{
		cout << v + 1 << "\t";
	}
	cout << endl;
	for(int v = 0; v < Max_Vtx; v++)
	{
		cout << score[v] << "\t";
	}
	cout << endl;
/*
	cout << "score: " << endl;
	for(int v = 0; v < Max_Vtx; v++)
	{
		if(vectex[v]) cout << score[v] << "\t\t";
	}
	cout << endl;
*/
}

int tabu( int Max_Iter )
{
     int i, j, k, l, bestlen = 0, am, am1,  ti, m1;
	long double ww, ww1, ww2;
     Iter = 0;
     clearGamma(); 
     while( 1 )
     {
        am = selectC0();
        if( am != -1 )
        {
            l = expand( am );
            Iter++;
            //if( Wbest == Waim )
			if(long_double_equals(Wbest, Waim))
               return Wbest;
        }
        else 
            break;
     }
#ifdef linux_mode
     times(&finish);
     double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
			double finish_time = timeGetTime() - start_time;
#endif
	 finish_time = round(finish_time * 100)/100.0;
//cout << "Wbest: " << Wbest << endl;
	if(long_double_greater(Wbest, lbest))
		{
			lbest=Wbest; len_W = len; real_solve1 = real_solve2; iter_best = M_iter + Iter;
/*
cout << "1. " << endl;
cout << "lbest: " << lbest << endl;
cout << "real_solve1: " << real_solve1 << endl;
*/
        for( i = 0; i < Max_Vtx; i++ )
        {
            TTbest[ i ] = Tbest[ i ];
        }
      

		}

	 if(finish_time>time_limit){M_iter+=Iter;
	 //if(Wbest>lbest)
				if(!verify())
				{
					cout << "there is something wrong" << endl;
					exit(1);
				}
	 //printf("%.2f	%d %d\n", real_solve1,lbest, M_iter);


			cout << "o " << LL(lbest) << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(M_iter) / 1000.0 / finish_time << endl;	

/*
for(int i = 0; i < Max_Vtx; i++)
{
//cout << TTbest[i] << endl;
	if(TTbest[i]) cout << i + 1 << endl;
}
*/
	 exit(0);
	}
     

     while( Iter < Max_Iter )
     {

        am = WselectC0();
        am1 = WselectC1();
        if( (am != -1) && (am1 != -1) )
        {
            ww = score[ C0[ am ] ];
            ww1 = score[ C1[ am1 ] ] + score[ BC[ C1[ am1 ] ] ];
        
            //if( ww > ww1 )
			if(long_double_greater(ww, ww1))
            {
                l = expand( am );
                
                Iter++;
                //if( Wbest == Waim )
				if(long_double_equals(Wbest, Waim))
                   return Wbest;
            }
            else
            {
                l = plateau( am1 );
                //if( Wbest == Waim )
				if(long_double_equals(Wbest, Waim))
                    return Wbest; 
                Iter++;
            }
        }
        else if( (am != -1) && (am1 == -1) )
        {
             l = expand( am );
             //if( Wbest == Waim )
			if(long_double_equals(Wbest, Waim))
               return Wbest;
                
             Iter++;
        }
        else if( (am == -1) && (am1 != -1) )
        {
             ti = Mumi_Weigt();
             m1 = cruset[ ti ];
             ww1 = score[ C1[ am1 ] ] + score[ BC[ C1[ am1 ] ] ];
             ww2 = score[ m1 ];
             //if( ww1 > ww2 )
			if(long_double_greater(ww1, ww2))
             {
                l = plateau( am1 );
                //if( Wbest == Waim )
				if(long_double_equals(Wbest, Waim))
                    return Wbest; 
                Iter++;
             }
             else
             {
                 k = backtract();
                 if( k == -1 )
                     return Wbest;
                 Iter++;
             }
        }
        else if( (am == -1) && (am1 == -1) )
        {
             k = backtract();
             if( k == -1 )
                return Wbest;
             Iter++;
        }
#ifdef linux_mode
        times(&finish);
        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
				double finish_time = timeGetTime() - start_time;
#endif
//show_state();
		if(long_double_greater(Wbest, lbest)){
					lbest=Wbest;
					real_solve1 = real_solve2;
					len_W = len;
/*
cout << "2. within a step:" << endl;
cout << "lbest: " << lbest << endl;
cout << "real_solve1: " << real_solve1 << endl;
*/

        for( i = 0; i < Max_Vtx; i++ )
        {
            TTbest[ i ] = Tbest[ i ];
        }

				}
		finish_time = round(finish_time * 100)/100.0;
        if(finish_time > time_limit){

			M_iter+=Iter;
        	//if(Wbest>lbest){
				if(!verify())
				{
					cout << "there is something wrong" << endl;
					exit(1);
				}
			//printf("%.2f	%d %d\n", real_solve1, lbest, M_iter);

			cout << "o " << LL(lbest) << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(M_iter) / 1000.0 / finish_time << endl;
/*
for(int i = 0; i < Max_Vtx; i++)
{
//cout << TTbest[i] << endl;
	if(TTbest[i]) cout << i + 1 << endl;
}
*/
			//exit(0);
			return Wbest;
			}
		}
     return Wbest;
}

inline int find_edge(int u, int v)
{
	for(int i = 0; i < neighbor_len[u]; i++)
	{
		int w = neighbor[u][i];
		int e = adj_edge[u][i];
		if(w == v) return e;
	}
	return -1;
}

int verify()
{
     int i, j, k1, k2, l, m;
		 long double verified_clique_weight = 0;
     for( i = 0; i < Max_Vtx; i++ )
     {
          if( TTbest[ i ] == 1 )
          {
//cout << i + 1 << " is in clique" << endl;
              for( j = i+1; j < Max_Vtx; j++ )
							{
		            if( TTbest[ j ] == 1 )
								{
									if(edge_is(i, j) == 1)
									{
				              cout << "hello there is something wrong" << endl;
											return 0;
									}
									int e = find_edge(i, j);
//cout << "the edge connecting " << i + 1 << " and " << j + 1 << " is " << e << endl;
									verified_clique_weight += We[e];
								}
							}
          }
     }
		if(!long_double_equals(lbest, verified_clique_weight))
		{
			cout << "lbest: " << lbest << endl;
			cout << "verified_clique_weight: " << verified_clique_weight << endl;
			cout << "weight calculated incorrectly" << endl;
			return 0;
		}
    // cout << "verified" << endl;
		return 1;
}

// WYY: Validate that the cruset is indeed a clique
int validate() {
	int i, j, k1, k2, l, m;
	long double verified_clique_weight = 0;
     for( i = 0; i < Max_Vtx; i++ )
     {
          if( vectex[ i ] == 1 )
          {
              for( j = i + 1; j < Max_Vtx; j++ )
							{
		            if( (vectex[ j ] == 1) && ( edge_is(i, j)== 1 ) ) {
		                cout << "hello there is something wrong" << endl;
		            	return 0;
		            }
								if(vectex[j] == 1)
								{
									int e = find_edge(i, j);
									verified_clique_weight += We[e];
//cout << "\t+= " << We[e] << endl;
								}
							}
          }
     }
		if(!long_double_equals(Wf, verified_clique_weight))
		{
			cout << "Wf: " << Wf << endl;
			cout << "verified_clique_weight: " << verified_clique_weight << endl;
			cout << "weight calculated incorrectly" << endl;
			return 0;
		}
     cout << "verified*****************************" << endl;		
	return 1;
}
/*
void Output()
{
    int i , j, k, l, sum; 
    FILE *fp ;
    int len = strlen(File_Name);
    
    strcpy(outfilename,File_Name) ;
    outfilename[len]='.';
    outfilename[len+1]='o';
    outfilename[len+2]='u';
    outfilename[len+3]='t';
    outfilename[len+4]='\0';

    fp = fopen(outfilename, "a+"); 
    for( i = 0; i < 1; i++ )
    {
        fprintf(fp, "sum = %6d, iter = %6d, len = %5d,  time = %8lf \n", 
        	W_used[ i ], Iteration[ i ], len_used[ i ],  time_used[ i ] ); 
    }
    
    fclose(fp); // WYY
    return;
    
    fprintf(fp, "\n\n the total information: \n");
    int wavg, iteravg, lenbb, success;
    wavg = iteravg = lenbb = success = 0;
    int best_v = 0;
    double timeavg = 0; 
    for( i = 0; i < 100; i++ )
		if( W_used[ i ] > best_v )
		{
			best_v = W_used[ i ];  
			lenbb = len_used[ i ];
		}
    
    int count = 0;
    fprintf(fp, "\n The best weight value for the maximum weighted problem is %6d \n", best_v);
    for( i = 0; i < 100; i++ )
    {
       wavg = wavg + W_used[ i ];
    }  
    double twavg = (double (wavg)) / 100 ; 
    for( i = 0; i < 100; i++ )
    if( W_used[ i ] == best_v )
    {
        count++;
        iteravg = iteravg + Iteration[ i ];
        timeavg = timeavg + time_used[ i ];
    }
    
    iteravg =  int ( (double (iteravg)) / count );
    timeavg = timeavg / (count*1000);
    fprintf(fp, "avg_sum = %10lf, succes = %6d, len = %5d, avg_iter = %6d,  time = %8lf \n", 
    			twavg, count, lenbb,  iteravg, timeavg );
    fclose(fp) ;
    return ;
}
*/
int Max_Tabu()
{
     int i, j, k, l, m;
     lbest = 0;
     int lenbest = 0;
     Titer = 0;
#ifdef linux_mode
	 times(&start);
#endif
#ifdef windows_mode
		start_time = timeGetTime();
#endif
	 while(1)
     {
//cout << 1 << endl;
         l = tabu(len_improve);
//cout << 2 << endl;
         M_iter = M_iter + Iter; 

         //if( l > lbest )
		if(long_double_greater(l, lbest))
         {
			 real_solve1=real_solve2;
			 lbest = l; 
			 len_W = len;
/*
cout << "3. after a whole try:" << endl;
cout << "lbest: " << lbest << endl;
cout << "real_solve1: " << real_solve1 << endl;
*/
			 Titer = M_iter; 
			 len_W = len_best;
//cout << "******update*****best****solution" << endl;  
        for( i = 0; i < Max_Vtx; i++ )
        {
            TTbest[ i ] = Tbest[ i ];
        }
      
         }
         
         //if( l >= Waim )
		if(long_double_geq(l, Waim))
           return lbest;
         //cout << " l = " << l << " i = " << i << endl;
		 
		 // wyy: clear configuration Information for the next restart
		 for(int j = 0; j < Max_Vtx; j++) {
          	conf_change[j] = 1;
          	time_stamp[j] = 0;
         }
#ifdef linux_mode
		 times(&finish);
		 double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
			unsigned long finish_time = timeGetTime() - start_time;
#endif
		 finish_time = round(finish_time * 100)/100.0;
		 if(finish_time>time_limit){
				if(!verify())
				{
					cout << "there is something wrong" << endl;
					exit(1);
				}
				//printf("%.2f	%d\n", real_solve1,lbest);
/*
			cout << "o " << LL(lbest) << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(M_iter) / 1000.0 / finish_time << endl;
*/
/*
for(int i = 0; i < Max_Vtx; i++)
{
//cout << TTbest[i] << endl;
	if(TTbest[i]) cout << i + 1 << '\t';
}
*/
			exit(0);
	}
     
     }
     return lbest;
}

	int get_best_weighted_clique_size()
	{
		return len_W;
	}

	void get_best_weighted_clique(int *outClique)
	{
		int i = 0;
		for(int j = 0; j < Max_Vtx; j++)
		{
			if(TTbest[j]) 
				//cout << j + 1 << endl;
				outClique[i++] = j + 1;
		}
	}

int main(int argc, char **argv)
{
     int seed;
	 File_Name = argv[1];
	 //Wmode=200;
	 time_limit=atof(argv[3]);
	 len_improve=4000;
	 seed = atoi(argv[2]);
	 //printf("%s	",argv[1]);
	 //Waim=INT_MAX;
	Waim = DBL_MAX;
	 srand(seed);
	 Initializing();
	 tm1 = Max_Vtx*sizeof( int );
     //  cout << "finish reading data" << endl;
     int i, l;
	 len_time = (int (100000000 / len_improve) ) + 1;
    //   cout << "len_time = " << len_time << endl;
     l = Max_Tabu();


}

#include <stdio.h>
#include <cassert>
#include "gtools.h"
#define MAXK 30
#define MAX_N 9
#define MAXNPATHS  10000

int session[MAXK][2];
int sessionMark[MAX_N];
int map[MAX_N][MAX_N];
int n,k;
int cutset[MAX_N];
int count = 0, nTarget = 0;
struct pathInfo{
	int n;
	int lens[MAXNPATHS];
	int a[MAXNPATHS][MAX_N];
}paths[MAXK];

void printNetwork()
{
	/* print the graph */
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++)
			printf("%d ", map[i][j]);
		printf("\n");
	}
	printf("with sessions:"); 
	for (int i=0; i<k; i++) printf("%d %d;", session[i][0], session[i][1]);
	printf("\n"); 
}
void reset(int *a, int k)
{
	for (int i=0; i<k; i++) a[i] = 0;
}
bool next(int *a, int k, int base) //return if there is a next one
{
	int i=0;
	for (i=0; i<k; i++) {
		if (a[i] < base-1) break;
		a[i] = 0;
	}
	if (i>=k) return false; //overflow
	a[i]++;
	return true;
}
void setSessions()
{
	k = 0;
	for (int i=0; i<n; i++)
		for (int j=i+1; j<n; j++) if (map[i][j]>1){
			session[k][0] = i;
			session[k][1] = j;
			k++;
		}
}
void enuPaths(int * path, bool * mark, int s, int t, int i, int m)
{
	int last = (i)? path[i-1] : s;
	if (last == t) {
		int idx = paths[m].n++;
		paths[m].lens[idx] = i;
		for (int j=0; j<i; j++) paths[m].a[idx][j] = path[j];
		return;
	}
	for (int j=0; j<n; j++) if (!mark[j] && map[last][j] == 1){
		path[i] = j;
		mark[j] = true;
		enuPaths(path, mark, s, t, i+1, m);
		mark[j] = false;
	}
	return ;
}
void setPaths()
{
	int path[MAX_N];
	bool mark[MAX_N];
	for (int i=0; i<n; i++) mark[i] = false;
	for (int m=0; m<k; m++) {
		int s = session[m][0];
		int t = session[m][1];
		paths[m].n = 0;
		mark[s] = true;
		enuPaths(path, mark, s, t, 0, m);
		mark[s] = false;
	}
	return ;
}
int getNcross(int m, int i, int * cutset)
{
	int ncross = 0;
	assert(i<paths[m].n);
	for (int j=0, last = session[m][0]; j<paths[m].lens[i]; j++) {
		if (cutset[last] != cutset[paths[m].a[i][j]]) ncross++;
		last = paths[m].a[i][j];
	}
	return ncross;
}
bool isMerger(int c1, int c2)
{
	int cut1[MAX_N], cut2[MAX_N];
	for (int i=0,testbit=1; i<n-1; i++, testbit<<=1) { 
		cut1[i] = (c1 & testbit)? 1:0; 
		cut2[i] = (c2 & testbit)? 1:0; 
	}
	cut1[n-1] = cut2[n-1] = 1; 
	//check if c1+c2 is compatible
	int count = 0;
	for (int m=0; m<k; m++) {
		int leastCrossShortest = n; //the longest paths has length n-1
		int leastCrossNon = n;
		for (int i=0; i<paths[m].n; i++) {
			int ncross = getNcross(m,i,cut1) + getNcross(m,i,cut2);
			if (paths[m].lens[i] == map[session[m][0]][session[m][1]]) {
				if (leastCrossShortest != n && leastCrossShortest != ncross) return false;
				leastCrossShortest = ncross;
			}
			else {
				if (leastCrossNon > ncross) leastCrossNon = ncross;
			}
		}
		if (leastCrossNon < leastCrossShortest) return false;
		count += leastCrossShortest>2;
	}
	return count<=1;
}
bool isOrth(int *cutset, int m)
{
	int s = session[m][0];
	int t = session[m][1];
	for (int i=0; i<paths[m].n; i++) if (paths[m].lens[i] == map[s][t]) {
		if (getNcross(m,i,cutset) > 1) return false;
	}
	return true;
}
bool no_orth_cut(sparsegraph * sg)
{
	reset(cutset, n-1);
	cutset[n-1] = 1;
	int idx = 0;
	do{
		bool orthToAll = true;
		int nSize = 0;
		for (int i=0; i<n; i++) nSize += cutset[i];
		if (nSize == n) continue;
		for (int i=0; i<k; i++) {
			if (!isOrth(cutset, i)) { orthToAll = false; }
		}
		idx++;
		if (orthToAll) {
			return false;
		}
	}while (next(cutset, n-1, 2));
	return true;
}
void second_check(sparsegraph * sg)
{
	int nCuts = (1<<(n-1))-1;
	for (int i=0; i<nCuts; i++) {
		for (int j=i+1; j<nCuts; j++) {
			if (isMerger(i,j)) {
				return;
			}
		}
	}
	nTarget++;
	printf("find one! graph %d (e=%d).\n", count, sg->nde/2);
	printNetwork();
}
void check_graph(sparsegraph * sg)
{
	n = sg->nv;
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			map[i][j] = 0;
	for (int i=0; i<n; i++) {
		for (int j=0; j<sg->d[i]; j++) {
			int ii = sg->e[sg->v[i]+j];
			map[i][ii] = 1;
		}
	}
	for (int k=0; k<n; k++)
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++) {
				if (i==j || k==j || i==k) continue;
				if (map[i][k] == 0 || map[k][j] == 0) continue;
				if (map[i][j] > 0 && map[i][j] <= map[i][k] + map[k][j]) continue;
				map[i][j] = map[i][k] + map[k][j];
			}
	setSessions();
	if (k<= 2) return; // as an orthgonal cut must exist for two sessions
	setPaths();
	if (no_orth_cut(sg)) {
		second_check(sg);
		//printNetwork();
	}
}
int main(int argc, char** argv)
{
	if (argc < 2) {
		printf("please indicate the file containing the list of graphs.\n");
		return 0;
	}
	FILE * fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("cannot open file %s.\n", argv[1]);
		return 0;
	}
	int nloops;
	sparsegraph * sg;
	while((sg = read_sg_loops(fp, NULL, &nloops)) != NULL) { 
		count++; 
		//printf("Checking graph %d (e=%d).\n", count, sg->nde/2);
		check_graph(sg);
	}
	printf("Checked %d graphs, found %d configurations.\n", count, nTarget);
	fclose(fp);
	return 0;
}

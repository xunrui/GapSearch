#include <stdio.h>
#include "gtools.h"
#define MAXK 8
#define MAX_N 9

int session[MAXK][2];
int sessionMark[MAX_N];
int map[MAX_N][MAX_N];
int n,k=3;
int cutset[MAX_N];
int count = 0, nTarget = 0;
bool oMark[1<<(MAX_N-1)][MAXK];

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
void initSessionEnu()
{
	for (int i=0; i<n; i++) sessionMark[i] = 0;
	session[0][0] = session[0][1] = -1;
}
bool nextSession(int j)
{
	int startS = session[j][0];
	int startT = session[j][1]+1;
	if (startS < 0) {
		startS = (j==0)? 0:session[j-1][0]+1;
		startT = startS + 1;
	}else {
		sessionMark[startS] = sessionMark[startT-1] = 0;
	}
	for (int u=startS; u<n; u++) if (sessionMark[u]==0) {
		int temp = (u==startS)? startT:u+1;
		for (int v = temp; v<n; v++) if (sessionMark[v]==0 && map[u][v]>1) {
			session[j][0] = u;
			session[j][1] = v;
			sessionMark[u] = sessionMark[v] = 1;
			if (j==k-1) return true;
			else {
				j++;
				session[j][0] = -1;
				return nextSession(j);
			}
		}
	}
	if (j==0) return false;
	return nextSession(j-1);
}

int cut1[MAX_N], cut2[MAX_N], path[MAX_N], currentBest, best;
bool searchPath(int i, int ncross, int m, sparsegraph * sg)
{
	int u = path[i-1];
	if (u == session[m][1]) {
		if (i-1 == map[session[m][0]][session[m][1]]) {
			if (best >= 0) {
				if (ncross != best) return false; //two shortest paths have different cross count.
			}
			else {
				if (currentBest > -1 && currentBest < ncross) return false; 
				// there is a non-shortest path has less cross count.
				best = ncross;
			}
		}
		else {
			if (best >= 0 && ncross < best) return false;
			if (currentBest > ncross || currentBest == -1) currentBest = ncross;
		}
		return true;
	}
	for (int iv = 0; iv<sg->d[u]; iv++) {
		int v = sg->e[sg->v[u]+iv];
		int delt=0;
		if (cut1[u] != cut1[v]) delt++;
		if (cut2[u] != cut2[v]) delt++;
		bool isInPath = false;
		for (int j=0; j<i-1; j++) 
			if (path[j] == v) { isInPath = true; break; }
		if (!isInPath) {
			path[i] = v;
			if (!searchPath(i+1, ncross+delt, m, sg)) return false;
		}
	}
	return true;
}
bool isMerger(int m, int c1, int c2, sparsegraph * sg)
{
	for (int i=0,testbit=1; i<n-1; i++, testbit<<=1) { 
		cut1[i] = (c1 & testbit)? 1:0; 
		cut2[i] = (c2 & testbit)? 1:0; 
	}
	cut1[n-1] = cut2[n-1] = 1; path[0] = session[m][0];
	currentBest = best = -1;
	if (searchPath(1, 0, m, sg)) {
		/*
		printf("s%d: ", m);
		for (int i=0; i<n; i++)
			if (cut1[i]) printf("%d ", i);
		printf("; ");
		for (int i=0; i<n; i++)
			if (cut2[i]) printf("%d ", i);
		printf("\n");
		*/
		return true;
	}
	return false;
}
bool isOrth(int m, sparsegraph * sg)
{
	int s = session[m][0];
	int t = session[m][1];
	int queue[MAX_N], crossTime[MAX_N], mark[MAX_N];
	int head=0, tail=1;
	for (int i=0; i<n; i++) { crossTime[i] = 0; mark[i] = 0;}
	queue[0] = t; mark[t] = 1;
	while (head<tail) {
		int u = queue[head];
		head++;
		for (int i=0; i<sg->d[u]; i++) {
			int v = sg->e[sg->v[u]+i];
			if (map[s][v] == map[s][u] - 1) {
				if (mark[v] == 0) {
					mark[v] = 1;
					queue[tail] = v;
					tail++;
				}
				int nCross = crossTime[u];
				if (cutset[u] ^ cutset[v]) {
					nCross++;
					if (nCross > 1) return false;
				}
				if (nCross > crossTime[v]) 
					crossTime[v] = nCross;
			}
		}
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
			oMark[idx][i] = isOrth(i,sg);
			if (!isOrth(i, sg)) { orthToAll = false; }
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
			bool merger = true;
			int sum = 0, m=-1;
			for (int o=0; o<k; o++)
				if (oMark[i][o] != oMark[j][o]) {
					merger = false;
					break;
				} else {
					sum += (oMark[i][o]) ? 1:0;
					if (!oMark[i][o]) m = o;
				}
			if (merger && sum == k-1) {
				//printf("%0o  %0o\n", i,j);
				if (isMerger(m,i,j,sg)) {
					//printf("Never mind. There is a merger.\n");
					return;
				}
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
	initSessionEnu();
	bool firstSearch = true;
	while (nextSession(firstSearch? 0:k-1)) {
		firstSearch = false;
		if (no_orth_cut(sg)) 
			second_check(sg);
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

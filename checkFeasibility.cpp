
void setCoef(int * coef, int m, int i)
{
	for (int j=0; j<nEdges; j++) coef[j] = 0;
	for (int j=0, last = session[m][0]; j<paths[m].lens[i]; j++) {
		int next = paths[m].a[i][j];
		coef[adj[last][next]] = 1;
		last = next;
	}
}
#define MAXNEDGES (MAX_N*MAX_N)
bool isFeasible()
{
	int nCons;
	int coef1[MAXNEDGES], coef2[MAXNEDGES];
	int ind[MAXNEDGES+1];
	double val[MAXNEDGES+1];


	glp_prob * lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_cols(lp, nEdges+1);
	for (int i=0; i<nEdges+1; i++) 
		glp_set_col_bnds(lp, i+1, GLP_LO, 0, 0);
	glp_set_obj_coef(lp, nEdges+1, 1);
	nCons = 0;
	for (int i=0; i<nEdges; i++) {
		nCons++;
		glp_add_rows(lp,1);
		glp_set_row_bnds(lp, nCons, GLP_UP, 0, 0);
		ind[1] = nEdges+1; val[1] = 1; //gamma <= l_e
		ind[2] = i+1; val[2] = -1;
		glp_set_mat_row(lp, nCons, 2, ind, val);
	}
	for (int m=0; m<k; m++) {
		int first = -1;
		for (int i=0; i<paths[m].n; i++) if (paths[m].isShort[i]) { first = i; break; }
		assert(first>=0);
		setCoef(coef1, m, first);
		for (int i=first+1; i<paths[m].n; i++) if (paths[m].isShort[i]) {
			setCoef(coef2, m, i);
			nCons++;
			glp_add_rows(lp,1);
			glp_set_row_bnds(lp, nCons, GLP_FX, 0, 0);
			int nEle = 0;
			for (int j=0; j<nEdges; j++) if (coef1[j] ^ coef2[j]) {
				nEle++;
				ind[nEle] = j+1;
				val[nEle] = (coef1[j])? 1:-1;
			}
			glp_set_mat_row(lp, nCons, nEle, ind, val);
		}
		for (int i=0; i<paths[m].n; i++) if (!paths[m].isShort[i]){
			setCoef(coef2, m, i);
			nCons++;
			glp_add_rows(lp,1);
			glp_set_row_bnds(lp, nCons, GLP_UP, 0, 0);
			int nEle = 0;
			for (int j=0; j<nEdges; j++) if (coef1[j] ^ coef2[j]) {
				nEle++;
				ind[nEle] = j+1;
				val[nEle] = (coef1[j])? 1:-1;
			}
			glp_set_mat_row(lp, nCons, nEle, ind, val);
		}
	}
	nCons++;
	glp_add_rows(lp,1);
	glp_set_row_bnds(lp, nCons, GLP_UP, 0, 1);
	for (int j=0; j<nEdges; j++) {
		ind[j+1] = j+1;
		val[j+1] = 1;
	}
	glp_set_mat_row(lp, nCons, nEdges, ind, val);
	glp_term_out(GLP_OFF);
	glp_simplex(lp, NULL);
	double ret = glp_get_obj_val(lp);
	glp_delete_prob(lp);
	return (ret>0);
}

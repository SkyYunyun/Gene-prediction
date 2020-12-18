#include<math.h>
#include"sBlast.h"

double SB_long_penalty = 0.02;
double SB_overlap_ratio = 0.45;
int SB_overlap_len = 1<<30; // that is overlap only be controlled by SB_overlap_ratio 
int SB_max_gap_len = 8192;
int SB_score_min_len = 11;

static double log_min_len;
static int *flag;

inline void VexEdge::add(node_t n)
{
	if (num == max) {
		max += MEM_BLOCK;
		array = (node_t*)realloc(array, sizeof(node_t) * max);
	}
	array[num++] = n;
}
inline bool SortBlast::isconsist(aBlock *n1, aBlock *n2) // a < b
{
	int tmp, tmp2;
	
	if (!(sb_flag & SB_INVERSION) && n1->dir != n2->dir) return false;
	tmp=n1->qp-n2->qt;
	if(tmp>0) {
		if(tmp>SB_overlap_len) return(false);
		if(double(tmp)/(n1->qp-n1->qt)>SB_overlap_ratio) return(false);
		if(double(tmp)/(n2->qp-n2->qt)>SB_overlap_ratio) return(false);
		tmp2 = 0;
	} else tmp2 = -tmp; // tmp2 is the half-gap length for the query
	if (sb_flag & SB_INVERSION) {
		if (sb_flag & SB_FORWARD) {
			if (n1->st > n2->st || n1->sp > n2->sp) return false;
			tmp = n1->sp - n2->st;
		} else {
			if (n2->st > n1->st || n2->sp > n1->sp) return false;
			tmp = n2->sp - n1->st;
		}
	} else {
		if (n1->dir) {
			if (n1->st > n2->st || n1->sp > n2->sp) return false;
			tmp = n1->sp - n2->st;
		} else {
			if (n2->st > n1->st || n2->sp > n1->sp) return false;
			tmp = n2->sp - n1->st;
		}
	}
	if(tmp>0) {
		if(tmp>SB_overlap_len) return(false);
		if(double(tmp)/(n1->sp-n1->st)>SB_overlap_ratio) return(false);
		if(double(tmp)/(n2->sp-n2->st)>SB_overlap_ratio) return(false);
		tmp = 0;
	} else tmp = -tmp; // tmp is the other half-gap length for subject
	if (!(sb_flag & SB_BOTH_GAP)) {
		if (tmp > SB_max_gap_len || tmp2 > SB_max_gap_len)
			return false;
	} else {
		if (tmp2 - tmp > SB_max_gap_len || tmp - tmp2 > SB_max_gap_len)
			return false;
	}
	return true;
}
double SortBlast::penalty(aBlock *n1, aBlock *n2) // n1->qt < n2->qt
{
	length_t ql, sl;
	
	ql = n2->qt - n1->qp;
	sl = (n1->st < n2->st) ? (n2->st - n1->sp) : (n1->st - n2->sp);
	ql = (ql > sl)? ql : sl;
	if (ql < 0) {
		return double((-ql)*2) / ((n1->score > n2->score)?
				  (n1->qp - n1->qt + n1->sp - n1->st)
				: (n2->qp - n2->qt + n2->sp - n2->st));
	}
	if (ql <= SB_score_min_len)
		return 0.0;
	return (log(ql) - log_min_len) * SB_long_penalty;
}
void SortBlast::build_list(void) // construct edge list for optimization
{
	node_t i, j, v, *q, *top, *r, tmp, *stack = new node_t[w_num];
	top=stack;
	for (i = 0; i < w_num; ++i)
		flag[i] = 0;
	for (i = 1; i < w_num; ++i) {
		tmp = wbs[i]->qp;
		for (j = i - 1; j >= 0; --j) {
			if (flag[j] != i && wbs[j]->qp < tmp && isconsist(wbs[j], wbs[i])) {
				ve[i].add(j);
				*top++=j;
				while(top>stack) { // find relevant vexes.
					v=*--top;
					flag[v]=i;
					for (q = ve[v].array; q < ve[v].num + ve[v].array; ++q)
						if (*(r = flag + *q) != i && *r < w_num) {
							*top++ = *q; *r = w_num;
						}
				}
			}
		}
	}
	for (i = 0; i < w_num; i++)
		if (!flag[i]) ve[w_num].add(i);
	delete[] stack;
}
score_t SortBlast::optimal(node_t *opt_ind) // find optimal path
{
	score_t opt, tmp;
	length_t *l;
	int j;
	node_t *q;
	VexEdge *p;
	
	l = new length_t[w_num + 1];
	for (j = 0; j <= w_num; ++j) {
		l[j] = 0; flag[j] = 0;
	}
	p = ve + w_num;
	for (q = p->array; q < p->num + p->array; ++q) {
		l[*q] = ve[*q].score;
		flag[*q] = w_num;
	}
	for (opt = 0, j = w_num - 1; j >= 0; --j) {
		p = ve + j;
		if (p->num == 0 && l[j] > opt) {
			opt = l[j];
			*opt_ind = j;
			continue;
		}
		for (q = p->array; q < p->array + p->num; ++q) {
			tmp = (ve[*q].score < p->score)? ve[*q].score : p->score;
			tmp = ve[*q].score + l[j] - score_t(tmp * penalty(wbs[*q], wbs[j]) + 0.5);
			if (tmp > l[*q]) { l[*q] = tmp; flag[*q] = j; }
		}
	}
	delete[] l;
	return opt;
}
int SortBlast::back_trace(node_t *stack,node_t opt_ind)
{
	node_t *top = stack;
	while (opt_ind < w_num) {
		*top++ = opt_ind;
		opt_ind=flag[opt_ind];
	}
	return(top-stack);
}
// compute the optimal path stored in stack whose length will be returned
// if the optimal path is on reverse strand, return value will be less than 0
node_t SortBlast::get_optimal(node_t *stack)
{
	node_t opt_ind, i, len;
	score_t sf;

	log_min_len = log(SB_score_min_len);
	flag = new node_t[w_num + 1]; // temp array
	for(i = 0; i < w_num; i++) {
		ve[i].rewind();
		ve[i].score = wbs[i]->score;
	}
	ve[w_num].score = 0;
	if (sb_flag & SB_INVERSION) {
		score_t sb;
		
		sb_flag |= SB_FORWARD; // set forward flag
		build_list(); // forward flag is useful in the function isconsist
		sf = optimal(&opt_ind);
		len = back_trace(stack, opt_ind);
		// if only one block is detected, the orientation should be
		// adjusted accordingly.
		if (len == 1 && !wbs[*stack]->dir) len = -len;
	
		for (i = 0; i < w_num; i++) ve[i].rewind();
		sb_flag &= ~SB_FORWARD; // set backward flag
		build_list(); // flag[] will be changed in this function
		sb = optimal(&opt_ind);
		if (sf < sb) { // backward result is larger, take it.
			sf = sb;
			len = -back_trace(stack, opt_ind); // if backward, len < 0
			if (len == -1 && wbs[stack[0]]->dir) len = -len;
		} // otherwise, take the forward one
	} else {
		build_list();
		sf = optimal(&opt_ind);
		len = back_trace(stack, opt_ind);
		if (!wbs[stack[0]]->dir) len = -len; // on reverse strand
	}
	// here sf will be the modified optimal score. but it is not returned.
	delete[] flag;
	return len;
}
score_t cal_part_score(aBlock *n1, aBlock *n2)
{
	length_t ql, sl;
	score_t sc;
	
	sc = n2->score;
	ql = (n1->qt < n2->qt)? (n1->qp - n2->qt + 1) : (n2->qp - n1->qt + 1);
	sl = (n1->st < n2->st)? (n1->sp - n2->st + 1) : (n2->sp - n1->st + 1);
	if (ql > sl) {
		if (ql > 0) { // exert overlap panelty
			if (n1->score < n2->score) n1 = n2;
			sc -= score_t(double(ql) / (n1->qp - n1->qt + 1) * n1->score + 0.5);
			return sc;
		}
	} else {
		if (sl > 0) {
			if (n1->score < n2->score) n1 = n2;
			sc -= score_t(double(sl) / (n1->sp - n1->st + 1) * n1->score + 0.5);
			return sc;
		}
	}
	return sc;
}
score_t SortBlast::cal_score(node_t *stack, node_t len)
{
	if (len == 0) return 0;
	node_t i;
	score_t whole;
	
	aBlock *n1, *n2;
	n1 = wbs[*stack];
	whole = n1->score;
	
	for (i = 1; i < len; i++) {
		n2 = wbs[stack[i]];
		whole += cal_part_score(n1, n2);
		n1 = n2;
	}
	return whole;
}

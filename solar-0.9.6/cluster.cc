#include "crossBlast.h"
#include "table2d.h"
#include "qsort.h"

int SB_max_cluster_gap = 100000000;
table2d<aBasicBlock, int> cluster_table;

typedef svector<aBlock*>::iterator aBlock_pt;
extern svector<aBlock*> wbs_backup;

static svector<int> sort_array;

struct SBCmpSubject
{
	inline bool operator() (const int a, const int b) const
	{ return wbs_backup[a]->st < wbs_backup[b]->st; };
};

static inline void sort_subject(int b_ind, int e_ind)
{
	svector<int>::iterator iter, iter2, start;
	int i;
	aBasicBlock bb, *p;
	
	sort_array.rewind();
	for (i = b_ind; i < e_ind; ++i)
		sort_array.push_back(i);
	quick_sort(sort_array.begin(), sort_array.size(), SBCmpSubject());
	
	p = wbs_backup[sort_array[0]];
	bb.st = p->st; bb.sp = p->sp;
	bb.qt = p->qt; bb.qp = p->qp; bb.score = p->score;
	start = sort_array.begin();
	for (iter = sort_array.begin(); iter < sort_array.end(); ++iter) {
		p = wbs_backup[*iter];
		if (p->st - bb.sp < SB_max_cluster_gap) {
			if (bb.sp < p->sp) bb.sp = p->sp;
			if (bb.qt > p->qt) bb.qt = p->qt;
			if (bb.qp < p->qp) bb.qp = p->qp;
			bb.score += p->score;
		} else {
			cluster_table.push_node(bb);
			for (iter2 = start; iter2 < iter; ++iter2)
				cluster_table.push_back(*iter2);
			bb.st = p->st; bb.sp = p->sp;
			bb.qt = p->qt; bb.qp = p->qp; bb.score = p->score; 
			start = iter;
		}
	}
	cluster_table.push_node(bb);
	for (iter2 = start; iter2 < iter; ++iter2)
		cluster_table.push_back(*iter2);
}
void sb_cluster()
{
	int start, qt, qp;
	register size_t i;
	aBasicBlock *p;

	cluster_table.rewind(); // clear the cluster_table
	qt = wbs_backup[0]->qt; qp = wbs_backup[0]->qp; // set original qt and qp
	start = 0; // set start as zero
	for (i = 0; i < wbs_backup.size(); ++i) {
		p = wbs_backup[i];
		if (p->qt - qp < SB_max_cluster_gap) {
			if (qp < p->qp) qp = p->qp; // update qp to longest one
		} else {
			sort_subject(start, i); // fill cluster_table
			qt = p->qt; qp = p->qp;
			start = i;
		}
	}
	sort_subject(start, i);
}

#include "sBlast.h"
#include "qsort.h"

int SB_rep_bound = 20;
int SB_rep_depth = 10;

struct SBTailStruct
{
	node_t node;
	length_t tail;
};
// Here, the confliction between "<" and ">" is not a bug.
inline bool operator< (const SBTailStruct &a, const SBTailStruct &b) { return a.tail > b.tail; }

// these lines are for comparison between aBlock. Note that there two kinds of rules (see also in sBlast.h)
struct SBCmpStart
{
	inline bool operator() (const aBlock &a, const aBlock &b) const { return a.qt < b.qt; };
};

SortBlast::SortBlast(FILE *f, unsigned t) : Blocks(f)
{
	ve_max = 0; ve = 0; sb_flag = t; w_num = 0;
}
SortBlast::~SortBlast(void)
{   
    for(int i = 0;i < ve_max; ++i) free(ve[i].array);
	    free(ve);
}

void SortBlast::reset() // virtual
{
	Blocks::reset();
	rep_regions.rewind();
	rbs.rewind();
	wbs.rewind();
}

void SortBlast::flag_repeat() // blocks must be sorted according to qt
{
	SBTailStruct *tails;
	node_t q_head, t_head;
	node_t tail_size;
	node_t i, j, tmp;
	double tmp_start, tmp_stop;
	length_t last_start, last_stop;
	svector<same_pair<length_t> >::iterator iter;
	aBlock *p;
	
	tails = new SBTailStruct[num];
	q_head = 0; last_start = last_stop = 0;
	quick_sort(blocks, num, SBCmpStart());
	
	// find repeat regions
	for (i = 1; i <= num; i++) {
		if (i != num && blocks[i].qt >= last_start && blocks[i].qp <= last_stop) {
			blocks[i].repeat = rep_regions.size() - 1; // flag repeat
		} else if (i == num || blocks[i].qt - blocks[q_head].qt >= SB_rep_bound) {
			if (i - q_head >= SB_rep_depth) {
				tail_size = i - q_head;
				for (j = 0; j < tail_size; j++) {
					tails[j].tail = blocks[j + q_head].qp;
					tails[j].node = j + q_head;
				}
				quick_sort(tails, tail_size);
				t_head = 0;
				for (j = 1; j < tail_size; j++) {
					if (tails[t_head].tail - tails[j].tail >= SB_rep_bound) {
						if (j - t_head >= SB_rep_depth) break;
							else while (tails[++t_head].tail - tails[j].tail >= SB_rep_bound);
					} 
				}
				if (j - t_head >= SB_rep_depth) { // a repeated region
					tmp_start = 0.0;
					tmp_stop = 0.0;
					for (j = t_head; j < tail_size; j++) {
						blocks[tails[j].node].repeat = rep_regions.size(); // flag
						tmp_start += blocks[tails[j].node].qt;
						tmp_stop += tails[j].tail;
					}
					last_start = length_t(tmp_start / (tail_size - t_head) + 0.5);
					last_stop = length_t(tmp_stop / (tail_size - t_head) + 0.5);
					if (rep_regions.size() > 0) {
						iter = rep_regions.begin() + rep_regions.size() - 1;
						if (last_start <= iter->y) {
							iter->y = last_stop; last_start = iter->x;
						} else rep_regions.push_back(same_pair<length_t>(last_start, last_stop));
					} else rep_regions.push_back(same_pair<length_t>(last_start, last_stop));
				}
			}
			if (i != num)
				while (blocks[i].qt - blocks[++q_head].qt >= SB_rep_bound || blocks[q_head].repeat >= 0);
		}
	}
	delete[] tails;

	// re-mask blocks
	tmp = SB_rep_bound / 2;
	if (!rep_regions.begin()) return;
	for (p = blocks, iter = rep_regions.begin(); p < blocks + num; ++p) {
		if (p->repeat >= 0) continue;
		while (p->qt >= iter->y + tmp) {
			++iter;
			if (iter == rep_regions.end()) return;
		}
		if (p->qt >= iter->x - tmp && p->qp <= iter->y + tmp)
			p->repeat = iter - rep_regions.begin(); // mask
	}
}

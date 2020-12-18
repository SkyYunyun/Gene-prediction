#include <stdio.h>
#include "crossBlast.h"
#include "qsort.h"
#include "table2d.h"

svector<aBlock*> wbs_backup;

static svector<aBlock*> fragments, frag_backup;
static int *sb_stack;
static int stack_size, max_s_size;

CrossBlast::CrossBlast(FILE *fp, unsigned t):SortBlast(fp, t) {}
CrossBlast::~CrossBlast(void) {}
// judge if node a is overlapped with b (a<b)
inline bool CrossBlast::isconflict(aBlock *n1, aBlock *n2)
{
	length_t tmp;

	if (!(sb_flag & SB_CROSS_SBJCT)) {
		tmp = n1->qp - n2->qt;
		if (tmp > 0) {
			if (tmp > SB_overlap_len) return true;
			if (double(tmp) / (n1->qp - n1->qt) > SB_overlap_ratio) return true;
			if (double(tmp) / (n2->qp - n2->qt) > SB_overlap_ratio) return true;
		}
	}
	tmp = (n1->st < n2->st)? (n1->sp - n2->st) : (n2->sp - n1->st);
	if (tmp > 0) {
		if (tmp > SB_overlap_len) return true;
		if(double(tmp) / (n1->sp - n1->st) > SB_overlap_ratio) return true;
		if(double(tmp) / (n2->sp - n2->st) > SB_overlap_ratio) return true;
	}
	return(false);
}
bool CrossBlast::gen_work_blocks(node_t *stack, node_t stack_size)
{
	int i, j;
	bool *flag = new bool[w_num];
	node_t *p;

	for (i = 0; i < w_num; i++)
		flag[i] = true;
	if (stack_size < 0)
		stack_size = -stack_size;
	for (p = stack; p < stack + stack_size; p++) {
		flag[*p] = false;
		for (i = 0; i < *p; i++)
			if (flag[i] && isconflict(wbs[i], wbs[*p]))
				flag[i] = false;
		for (i = *p + 1; i < w_num; i++)
			if (flag[i] && isconflict(wbs[*p], wbs[i]))
				flag[i] = false;
	}
	for (j = 0, i = 0; i < w_num; i++) {
		if (flag[i])
			wbs[j++] = wbs[i];
	}
	w_num = j;
	wbs.resize(w_num);
	delete[] flag;
	if (w_num) return true;
		return false;
}
// generate all the working blocks for the same query and same subject
int CrossBlast::gen_first_work_blocks(int *c_site)
{
	register int i;
	int c_sbjct;
	
	if (*c_site == num) return -1;
	c_sbjct = blocks[*c_site].s_locus_id;
	wbs.rewind();
	if (sb_flag & SB_NOT_MASK_REPEAT) { // do not mask repeat
		for (i = *c_site; i < num && blocks[i].s_locus_id == c_sbjct; ++i)
			wbs.push_back(blocks + i);
	} else { // mask repeat
		for (i = *c_site; i < num && blocks[i].s_locus_id == c_sbjct; ++i) {
			if (blocks[i].repeat < 0) {
				wbs.push_back(blocks + i);
			} else rbs.push_back(blocks + i);
		}
	}
	w_num = wbs.size();
	*c_site = i;

	// allocate enough memory for ve
	if (w_num + 1 >= ve_max) {
		i = ve_max;
		ve_max = w_num + 1;
		ve = (VexEdge*)realloc(ve, sizeof(VexEdge) * ve_max);
		// realloc do not execute construction-function automatically,
		// so I have to do it manually.
		for (; i < ve_max; i++) {
			ve[i].array = 0;
			ve[i].num = ve[i].max = 0;
		}
	}
	return c_sbjct;
}
static void handle_fragment(FILE *fpout)
{
	svector<aBlock*>::iterator iter;
	fragments.rewind();
	for (iter = frag_backup.begin(); iter < frag_backup.end(); ++iter)
		if ((*iter)->repeat == SB_UNALIGNED)
			fragments.push_back(*iter);
	if (fragments.size() == 0) { // no fragments
		fprintf(fpout, "F\t0\n");
	} else {
		fprintf(fpout, "F\t%d\t", fragments.size());
		for (iter = fragments.begin(); iter < fragments.end(); ++iter)
			fprintf(fpout, "%d,%d;", (*iter)->qt, (*iter)->qp);
		fputc('\t', fpout);
		for (iter = fragments.begin(); iter < fragments.end(); ++iter) {
			if ((*iter)->dir)
				fprintf(fpout, "%d,%d;", (*iter)->st, (*iter)->sp);
			else
				fprintf(fpout, "%d,%d;", (*iter)->sp, (*iter)->st);
		}
		fputc('\t', fpout);
		for (iter = fragments.begin(); iter < fragments.end(); ++iter)
			fprintf(fpout, "%c%d;", ((*iter)->dir) ? '+' : '-', (*iter)->score);
		fputc('\n', fpout);
	}
}
// note that svector<int>::iterator is just int*
void CrossBlast::gen_alignment(FILE *fpout, SBResultStruct *rs, int sbjct)
{
	if (!(sb_flag & SB_BRIEF) && w_num == 1) { // just one block, output
		wbs[0]->repeat = SB_ALIGNED;
		if (wbs[0]->dir)
			fprintf(fpout, "A\t%d\t%d\t+\t%d\t%d\t1\t%d\t%d,%d;\t%d,%d;\t+%d;\n",
					wbs[0]->qt, wbs[0]->qp, wbs[0]->st, wbs[0]->sp, wbs[0]->score,
					wbs[0]->qt, wbs[0]->qp, wbs[0]->st, wbs[0]->sp, wbs[0]->score);
		else
			fprintf(fpout, "A\t%d\t%d\t-\t%d\t%d\t1\t%d\t%d,%d;\t%d,%d;\t-%d;\n",
					wbs[0]->qt, wbs[0]->qp, wbs[0]->st, wbs[0]->sp, wbs[0]->score,
					wbs[0]->qt, wbs[0]->qp, wbs[0]->sp, wbs[0]->st, wbs[0]->score);
		if (sb_flag & SB_MULTI) {
			fprintf(fpout, "T\t1\t%d\n", wbs[0]->score);
			fprintf(fpout, "F\t0\n");
		}
	} else { // two or more blocks will trigger dynamic programing.
		int score = 0, count = 0;
		do {
			if (w_num > max_s_size) {
				max_s_size = w_num;
				sb_stack = (int*)realloc(sb_stack, sizeof(int) * max_s_size);
			}
			stack_size = get_optimal(sb_stack);
			if (!write_result(sb_stack, stack_size, rs)) break; // too short
			score += rs->score; ++count;
			output(*rs, fpout, sbjct);
		} while ((sb_flag & SB_MULTI) && gen_work_blocks(sb_stack, stack_size));
		if ((sb_flag & SB_MULTI) && !(sb_flag & SB_BRIEF)) {
			fprintf(fpout, "T\t%d\t%d\n", count, score);
			handle_fragment(fpout); // output fragments
		}
	}
}
void CrossBlast::handle(FILE *fpout) // similar to handle()
{
	extern table2d<aBasicBlock, int> cluster_table;
	aBasicBlock *r;
	int i, c_site, c_sbjct;
	SBResultStruct rs;
	svector<int> *q;
	svector<int>::iterator iter;
	
	stack_size = 0; max_s_size = 0; sb_stack = 0;
	rs.list = 0; rs.max = 0;
	while (read_blocks()) {
		if (!(sb_flag & SB_NOT_MASK_REPEAT))
			flag_repeat();
		if (!(sb_flag & SB_BRIEF)) {
			fprintf(fpout, "Q\t%s\t%d\n", q_locus, q_len);
			if (rep_regions.size()) {
				fprintf(fpout, "R\t%d\t", rep_regions.size());
				for (size_t i = 0; i < rep_regions.size(); i++)
					fprintf(fpout, "%d,%d;", rep_regions[i].x, rep_regions[i].y);
				fprintf(fpout, "\n");
			} else fprintf(fpout, "R\t0\n");
		}
		if ((sb_flag & SB_NOT_MASK_REPEAT) || s_locus.size() > 1)
		// Otherwise, there is no need to perform sort another time
			quick_sort(blocks, num);

		c_site = 0;
		while ((c_sbjct = gen_first_work_blocks(&c_site)) >= 0) {
			if (!w_num) continue; // fully masked
			if (!(sb_flag & SB_BRIEF)) fprintf(fpout, "S\t%s\t%d\n", s_locus[c_sbjct], s_len[c_sbjct]);
			wbs_backup.rewind();
			for (i = 0; i < w_num; ++i)
				wbs_backup.push_back(wbs[i]);
			if (sb_flag & SB_CLUSTER) {
				sb_cluster(); // heirarchical clustering
				for (i = 0; i < int(cluster_table.size()); ++i) {
					q = cluster_table(i);
					r = &cluster_table[i];
					if (r->score < SB_min_chain_score) continue;
					// output cluster information
					if (!(sb_flag & SB_BRIEF))
						fprintf(fpout, "C%.6d\t%d\t%d\t*\t%d\t%d\t%d\t%d\n", i + 1,
								r->qt, r->qp, r->st, r->sp, q->size(), r->score);
					// generate working blocks
					wbs.rewind(); frag_backup.rewind();
					q->sort();
					for (iter = q->begin(); iter < q->end(); ++iter) {
						wbs.push_back(wbs_backup[*iter]);
						frag_backup.push_back(wbs_backup[*iter]);
					}
					w_num = wbs.size();
					// output alignment by dynamic programming
					gen_alignment(fpout, &rs, c_sbjct);
				}
			} else gen_alignment(fpout, &rs, c_sbjct);
		}
	}
	free(sb_stack);
	free(rs.list); // I have to free this manually.
}

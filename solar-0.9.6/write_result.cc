#include <stdio.h>
#include "sBlast.h"
#include "qsort.h"

int SB_min_chain_score = 0;

bool SortBlast::write_result(node_t *stack, node_t l, SBResultStruct *rs)
{
	node_t i;

	if (l < 0) { rs->dir = false; rs->num = -l; }
		else { rs->dir = true; rs->num = l; }
	rs->score = cal_score(stack, rs->num);
	if (rs->score < SB_min_chain_score)
		return false;
	rs->qt = wbs[*stack]->qt;
	rs->qp = wbs[stack[rs->num - 1]]->qp;
	if (l >= 0) {
		rs->st = wbs[*stack]->st;
		rs->sp = wbs[stack[rs->num - 1]]->sp;
	} else {
		rs->st = wbs[stack[rs->num - 1]]->st;
		rs->sp = wbs[*stack]->sp;
	}
	if (rs->num > rs->max) {
		rs->max = rs->num;
		rs->list = (aBlock**)realloc(rs->list, sizeof(aBlock*) * rs->max);
	}
	for (i = 0; i < rs->num; i++) {
		rs->list[i] = wbs[stack[i]];
		wbs[stack[i]]->repeat = SB_ALIGNED;
	}
	return true;
}
void SortBlast::output(const SBResultStruct &rs, FILE *fpout, int sbjct)
{
	node_t i;
	aBlock *p;
	if (sb_flag & SB_BRIEF) {
		fprintf(fpout, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t", q_locus, q_len, rs.qt, rs.qp, (rs.dir)?'+':'-',
				s_locus[sbjct], s_len[sbjct]);
	} else fprintf(fpout, "A\t%d\t%d\t%c\t", rs.qt, rs.qp, (rs.dir)?'+':'-');
	fprintf(fpout, "%d\t%d\t%d\t%d\t", rs.st, rs.sp, rs.num, rs.score);
	
	for (i = 0; i < rs.num; i++)
		fprintf(fpout, "%d,%d;", rs.list[i]->qt, rs.list[i]->qp);
	fprintf(fpout, "\t");
	for (i = 0; i < rs.num; i++) {
		p = rs.list[i];
		if (p->dir) // orientation
			fprintf(fpout, "%d,%d;", p->st, p->sp);
		else fprintf(fpout, "%d,%d;", p->sp, p->st);
	}
	fprintf(fpout, "\t");
	for (i = 0; i < rs.num; i++) {
		if (rs.list[i]->dir) fputc('+', fpout);
			else fputc('-', fpout);
		fprintf(fpout, "%d;", rs.list[i]->score);
	}
	fputc('\n', fpout);
}

#include<string.h>
#include<stdlib.h>
#include"sBlast.h"
#include"qsort.h"

unsigned sb_flag = 0;

Blocks::Blocks(FILE *f)
{
	fp=f;num=max=0;blocks=0;
	curr=(aLine*)malloc(sizeof(aLine));
	last=(aLine*)malloc(sizeof(aLine));
	q_locus=(char*)malloc(LOCUS_LEN);
	read_line(curr); // read the first record
	max_send = 0;
}
Blocks::~Blocks(void)
{
	free(curr);free(last);free(blocks);
	free(q_locus);
	Blocks::reset();
}
// reset all the variable for the next query
void Blocks::reset()
{
	num = 0; // clear all the blocks
	for (size_t i = 0; i < s_locus.size(); ++i)
		delete[] s_locus[i]; // free the memory for all the subject locus
	s_locus.rewind();
	s_len.rewind(); // s_len must be synchronized with s_locus anywhere
	q_len = 0;
}
// read a line from data file
void Blocks::read_line(aLine *line)
{
	char text[LOCUS_LEN], *q;
	int c = 0, i;

	line->num[0]=1;
	for(i = 0; i < 10; ++i) {
		if (i == 0) q = line->locus[0];
			else if (i == 2) q=line->locus[1];
				else q = text;
		while (!feof(fp) && (c = fgetc(fp)) != '\t' && c != ' ' && c != '\n')
			*q++ = c;
		if(feof(fp)) { line->num[0]=0;return; } // read final line.
		*q='\0';
		if (i == 9) {
			line->e = atof(text);
		} else if (i != 0 && i != 2) line->num[i] = atoi(text);
	}
	while (c != '\n') c = fgetc(fp);
}
void Blocks::add_first() // add the first blocks of a query
{
	aLine *tmp;
	char *tmp2;

	reset(); // clear all the previous records
	// set query name and length
	strcpy(q_locus, curr->locus[0]);
	q_len = curr->num[1];
	// set first subject name and length
	tmp2 = new char[strlen(curr->locus[1]) + 1]; // this will be freeed in reset()
	strcpy(tmp2, curr->locus[1]);
	s_locus.push_back(tmp2);
	s_len.push_back(curr->num[3]);
	add(curr->num, curr->e); // add the first block
	tmp = last; last = curr; curr = tmp; // change curr as last
}
node_t Blocks::read_blocks()
{
	if (curr->num[0] == 0) { num = 0; return 0; } // have read all the records in the file
	
	aLine *tmp;
	char *tmp2;

	add_first(); // add the current record that has not been added
	while (1) {
		read_line(curr);
		if (curr->num[0] == 0) break; // end of file
		if (strcmp(curr->locus[0], last->locus[0])) break; // curr is a new query
		if (strcmp(curr->locus[1], last->locus[1])) { // curr is a new subject
			if (last->num[3] == 0) s_len[s_len.size() - 1] = max_send;
			max_send = 0;
			tmp2 = new char[strlen(curr->locus[1]) + 1];
			strcpy(tmp2, curr->locus[1]);
			s_locus.push_back(tmp2);
			s_len.push_back(curr->num[3]);
		}
		add(curr->num, curr->e);
		tmp = curr; curr = last; last = tmp;
	}
	// subject has been changed due to the change of query, check length
	if (last->num[3] == 0) s_len[s_len.size() - 1] = max_send;
	max_send = 0;
	return num;
}
// add a block
void Blocks::add(int n[9], double e)
{
	if(num==max) {
		max+=MEM_BLOCK;
		blocks=(aBlock*)realloc(blocks,sizeof(aBlock)*max);
	}
	aBlock *tmp=blocks+num;
//	tmp->e = e;
	tmp->qt=n[4];tmp->qp=n[5];
	tmp->st=n[6];tmp->sp=n[7];
	tmp->score=n[8];
	tmp->repeat = SB_UNALIGNED; // SB_UNALIGNED (<0) stands for non-repeat
	
	int tmp2; // set the strand
	if (tmp->qt > tmp->qp) {
		tmp2 = tmp->qt; tmp->qt = tmp->qp; tmp->qp = tmp2;
		if (tmp->st > tmp->sp) {
			tmp->dir = true;
			tmp2 = tmp->st; tmp->st = tmp->sp; tmp->sp = tmp2;
		} else tmp->dir = false;
	} else {
		if (tmp->st > tmp->sp) {
			tmp->dir = false;
			tmp2 = tmp->st; tmp->st = tmp->sp; tmp->sp = tmp2;
		} else tmp->dir = true;
	}

	// adjust the estimation of query and subject length
	if (n[5] > q_len) q_len = n[5];
	if (tmp->sp > max_send) max_send = tmp->sp;
	
	tmp->s_locus_id = s_locus.size() - 1;
	++num;
}

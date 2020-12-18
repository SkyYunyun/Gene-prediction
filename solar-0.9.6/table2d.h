#ifndef TABLE2D_H_
#define TABLE2D_H_

#include "svector.h"

template <class TYPEN, class TYPEV>
class table2d 
{
public:
	typedef svector<TYPEV> vector_t;
private:
	size_t num, max;
	TYPEN *content;
	vector_t *list;
	vector_t *curr;
public:
	table2d() { construct(); };
	~table2d() { destroy(); };
	void construct() { content = 0; list = 0; num = 0; max = 0; curr = 0; };
	void destroy()
	{
		for (size_t i = 0; i < max; ++i)
			list[i].destroy();
		free(content); free(list);
	};
	inline void push_back(const TYPEV &a) { curr->push_back(a); };
	inline void push_back(size_t node, const TYPEV &a) { list[node].push_back(a); };
	inline TYPEV *push_null() { return curr->push_null(); };
	inline TYPEV *push_null(size_t node) { return list[node].push_null(); };
	inline void push_node(const TYPEN &a)
	{
		if (num == max) {
			max += LIH_BLOCK_SIZE;
			list = (vector_t*)realloc(list, sizeof(vector_t) * max);
			content = (TYPEN*)realloc(content, sizeof(TYPEN) * max);
			// realloc do not run constructor, so I have to run it by myself.
			for (size_t i = num; i < max; ++i)
				list[i].construct();
		}
		content[num] = a;
		curr = list + num;
		++num;
	}
	inline vector_t *push_node(TYPEN **top)
	{
		if (num == max) {
			max += LIH_BLOCK_SIZE;
			list = (vector_t*)realloc(list, sizeof(vector_t) * max);
			content = (TYPEN*)realloc(content, sizeof(TYPEN) * max);
			// realloc do not run constructor, so I have to run it by myself.
			for (size_t i = num; i < max; ++i)
				list[i].construct();
		}
		curr = list + num;
		*top = content + num;
		++num;
		return curr;
	}
	inline void set(size_t i) { curr = list + i; };
	inline TYPEN &operator[] (size_t i) { return content[i]; };
	inline vector_t *operator() (size_t i) { return list + i; };
	inline size_t size() { return num; };
	inline size_t capacity() { return max; };
	inline void rewind()
	{
		for (size_t i = 0; i < num; ++i)
			list[i].rewind();
		num = 0;
	}
	inline void resize(size_t new_size)
	{
		if (new_size < max) {
			num = new_size;
		} else {
			list = (vector_t*)realloc(list, sizeof(vector_t) * new_size);
			content = (TYPEN*)realloc(list, sizeof(TYPEN) * new_size);
			// realloc do not run constructor, so I have to run it by myself.
			for (size_t i = max; i < new_size; ++i)
				list[i].construct();
			max = new_size; num = new_size;
		}
	}
};

#endif

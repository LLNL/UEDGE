

#ifndef LINKAGE

#define BLOCK_ALL 0xfffff

struct linkage {
	void * previous;
	void * next;
};
#define LINKAGE 1
#endif

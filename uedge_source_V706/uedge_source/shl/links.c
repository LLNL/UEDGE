
/***************************************************/
/*      Copyright (c) 1988-1994                    */
/*      Regents of the University of California    */
/*      Lawrence Livermore National Laboratory     */
/*      All rights reserved                        */
/***************************************************/

static unsigned int smask;
#include <signal.h>
#include "links.h"

struct linkage *add_link(item,head,tail)
struct linkage *item;
struct linkage **head;
struct linkage **tail;
{
    struct linkage *temp;
    temp = *tail;
#ifndef vms
		smask = sigsetmask(BLOCK_ALL);
#endif
		if(item == 0){
    } else if(*head == 0){
        *head=item;
        item->previous=0;
        item->next=0;
    } else {
        temp->next = item;
        item->previous = *tail;
        item->next=0;
    }
    *tail=item;
#ifndef vms
		sigsetmask(smask);
#endif
    return(item);
}

void rem_link(item,head,tail)
struct linkage *item;
struct linkage **head;
struct linkage **tail;
{
    struct linkage *before,*after;
#ifndef vms
		smask = sigsetmask(BLOCK_ALL);
#endif
    before = item->previous;
    after = item->next;
    if(before == 0 && after == 0){
        *head = 0;
        *tail = 0;
    } else if(item == *head){
        *head = after;
        before = *head;
        before->previous = 0;
    } else if (item == *tail){
        *tail=item->previous;
        after = *tail;
        after->next = 0;
    } else {
        before = item->previous;
        after = item->next;
        before->next=after;
        after->previous = before;
    }
#ifndef vms
		sigsetmask(smask);
#endif
}


struct linkage *insert_link(item,ins_item,head,tail)
struct linkage *item,*ins_item;
struct linkage **head;
struct linkage **tail;
{
    struct linkage *s_next,*s_prev;
#ifndef vms
		smask = sigsetmask(BLOCK_ALL);
#endif
    s_next = item->next;
    s_prev = item->previous;
		if(ins_item == 0){
		} else if(*head == 0){
			add_link(ins_item,head,tail);
		} else if (item == 0){
			add_link(ins_item,head,tail);
    } else if(*head == item){
        *head = ins_item;
        ins_item->previous = 0;
        ins_item->next = item;
        item->previous = ins_item;
    } else {
        s_prev->next = ins_item;
        ins_item->previous = s_prev;
        ins_item->next = item;
        item->previous = ins_item;
    }
#ifndef vms
		sigsetmask(smask);
#endif
		return(ins_item);
}

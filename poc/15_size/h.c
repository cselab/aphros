#include <search.h>
#include <stdio.h>
#include "h.h"

enum { N = 1024 };
static char h_str[N];
int
hash(float key, char *s)
{
    sprintf(s, "%ld", (size_t)key);
    return 0;
}

int
h_ini(int n)
{
    return hcreate(n);
}

int
h_fin(void)
{
    hdestroy();
    return 0;
}

int
h_find(float key)
{
    ENTRY e, *p;
    hash(key, h_str);
    e.key = h_str;
    p = hsearch(e, FIND);
    if (p == NULL) {
	return -1;
    } else
	return p->data - NULL;
}

void
h_enter(float key, int val)
{
    ENTRY e;
    
    hash(key, h_str);
    e.key = h_str;
    e.data = NULL + val;
    hsearch(e, ENTER);
}

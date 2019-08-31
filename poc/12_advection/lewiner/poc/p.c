#include <stdio.h>
#define C(i) (i ? 'I' : 'O')
int main() {
    int i, j, k, p;
    for (p = 0; p < 8; ++p) {
	i = ((p ^ (p >> 1)) & 1);
	j = ((p >> 1) & 1);
	k = ((p >> 2) & 1);
	printf("%c%c%c,\n", C(i), C(j), C(k));
    }
}

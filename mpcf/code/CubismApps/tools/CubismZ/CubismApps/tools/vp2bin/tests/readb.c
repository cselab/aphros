#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	FILE *fp;
	float fl;

	if (argc != 2) {
		printf("usage: %s filename\n", argv[0]);
		exit(1);
	}
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("file not found!\n");
		exit(1);
	}

	while (fread(&fl, 1, 4, fp) == 4) {
	        /* byte swap here */
		printf("%f\n", fl);
	}

	fclose(fp);
}

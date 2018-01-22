make:
	gcc src/zopfli/*.c -O2 -W -Wall -Wextra  -lm -o zopfli
#	gcc src/zopfli/*.c -O2 -W -Wall -Wextra -ansi -pedantic -lm -o zopfli

debug:
	gcc src/zopfli/*.c -g3 -lm -o zopfli

lib:
	gcc -c src/zopfli_lib/*.c -O2 -W -Wall -Wextra -ansi -pedantic
#	gcc -c src/zopfli_lib/*.c -g -ggdb  -W -Wall -Wextra -ansi -pedantic
#	ar rs libzopfli.a src/zopfli_lib/*.o
	ar rs libzopfli.a *.o
	ranlib libzopfli.a
	ln -sf src/zopfli_lib/zopfli.h .
	rm -f *.o

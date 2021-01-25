#!/bin/bash

prefix=_bloat_test_tmp_
numTranslationUnits=100

template='
#include <aphros/util/format.h>
#include <iostream>

void doFormat_a()
{
    const char* str1 = "somefile.cpp";
    const char* str2 = "asdf";
    std::cout << util::Format("{}", str1);
    std::cout << util::Format("{}{}", str1, 42);
    std::cout << util::Format("{}{}{}", str1, 42, str2);
    std::cout << util::Format("{}{}{}{}", str1, 42, 1, str2);
    std::cout << util::Format("{}{}{}{}{}", str1, 42, 1, 2, str2);
}
'

# Generate all the files
out=main.cpp
> $out
echo -e "#include \"all.h\"" >> $out
echo -e "#include <iostream>" >> $out
echo '
int main()
{' >> $out

for ((i=0;i<$numTranslationUnits;i++)) ; do
    n=$(printf "%03d" $i)
    f=${prefix}$n.cpp
    echo "$template" | sed -e "s/doFormat_a/doFormat_a$n/" -e "s/42/$i/" > $f
    echo "doFormat_a$n();" >> $out
    echo "void doFormat_a$n();" >> all.h
done

echo "return 0; }" >> $out

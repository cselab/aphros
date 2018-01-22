touch tags.lst
find | grep "\.c$" >> tags.lst
find | grep "\.cpp$" >> tags.lst
find | grep "\.h$" >> tags.lst
cscope -i tags.lst

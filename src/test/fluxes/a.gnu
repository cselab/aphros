set terminal x11 noenhanced
ll=system('echo ref_*')
pp="7 5 2 4 3"
p for [i=1:words(ll)] word(ll,i).'/p' w lp pt int(word(pp,i)) ps 2 lw 3 lt i t word(ll,i) , for [i=1:words(ll)] word(ll,i).'/vx' w lp pt int(word(pp,i)) ps 2 lw 3 lt i t word(ll,i) axis x1y2

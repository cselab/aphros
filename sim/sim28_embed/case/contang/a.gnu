set terminal qt noenhanced
plot for [f in system("ls stat*.dat")] f u "t":"m2" w l lw 2 t f

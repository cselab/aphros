include(lib.js)
changequote()

var nx = 5, ny = 5, u = []
InitGrid(5, 5, u)

function expo(x) {
  return Number.parseFloat(x).toExponential(16);
}

function print(x) { process.stdout.write(x) }
for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
	if (i > 0) print(" ")
	print(expo(u[i][j]))
    }
    print("\n")
}


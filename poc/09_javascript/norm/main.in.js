include(lib.js)
changequote()

var nx = 5, ny = 5, u = []
InitGrid(5, 5, u)

function expo(x) {
  return Number.parseFloat(x).toExponential(1);
}


function print(x) { process.stdout.write(x) }

var px = [], py = [], n = [];
Norm(2, 2, u, /**/ n)

for (j = -1; j < ny; j++) {
    for (i = -1; i < nx; i++) {
    print(expo(u[i][j]))
    print(" ")
    }
    print("\n")
}

print("========\n")
print(expo(n[0]))
print(" ")
print(expo(n[1]))
print("\n")

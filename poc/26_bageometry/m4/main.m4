define(`foreach_dimension',
`dnl
$1dnl
define(`x', `y')$1undefine(`x')dnl
define(`x', `z')$1undefine(`x')dnl
')dnl
foreach_dimension(
    n1.x = n.x*(b.x - a.x);
)dnl

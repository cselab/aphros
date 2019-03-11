changequote()

const fs = require('fs')
const readline = require('readline')

const rl = readline.createInterface({
    input: process.stdin
})

rl.on('line', (line) => {
    process.stdout.write(line)
    process.stdout.write("\n")
})

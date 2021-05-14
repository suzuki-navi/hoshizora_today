const fs = require('fs');

const sourceLines = fs.readFileSync("/dev/stdin", 'utf8').split("\n");

const startTime = process.argv[2];
const endTime = process.argv[3];

var is = 0;
var ip = 0;
var removeTime = "";
var removeBodyList = [];
while (true) {
    var sourceLine = sourceLines[is];
    if (sourceLine == "") {
        break;
    }
    const [sourceTime, sourceMessage] = sourceLine.split(" ", 2);
    if (sourceTime >= startTime && sourceTime < endTime) {
        console.log(sourceLine);
    }
    is++;
}



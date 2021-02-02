const fs = require('fs');

const sourceLines = fs.readFileSync("/dev/stdin", 'utf8').split("\n");
const patchLines = fs.readFileSync("../var/diff.txt", 'utf8').split("\n");

const startTime = process.argv[2];
const endTime = process.argv[3];

var is = 0;
var ip = 0;
var removeTime = "";
var removeBodyList = [];
while (true) {
    var patchLine = patchLines[ip];
    var sourceLine = sourceLines[is];
    if (sourceLine == "") {
        break;
    }
    if (patchLine == "") {
        patchLine = "9999-99-99T99:99 +";
    }
    const [patchTime, patchBody] = patchLine.split(" ", 2);
    const [sourceTime, sourceMessage] = sourceLine.split(" ", 2);
    if (sourceTime < patchTime) {
        if (sourceTime == removeTime) {
            var f = false;
            removeBodyList2 = [];
            for (const el of removeBodyList) {
                if (el == sourceMessage) {
                    f = true;
                } else {
                    removeBodyList2.push(el);
                }
            }
            if (!f) {
                if (sourceTime >= startTime && sourceTime < endTime) {
                    console.log(sourceLine);
                }
            }
            removeBodyList = removeBodyList2;
            if (removeBodyList.length == 0) {
                removeTime = "";
            }
        } else {
            if (removeTime != "") {
                for (const el of removeBodyList) {
                    if (removeTime >= startTime && removeTime < endTime) {
                        console.log(removeTime + " #-" + el);
                    }
                }
                removeTime = "";
                removeBodyList = [];
            }
            if (sourceTime >= startTime && sourceTime < endTime) {
                console.log(sourceLine);
            }
        }
        is++;
        continue;
    }
    if (patchBody.startsWith("+")) {
        if (patchTime >= startTime && patchTime < endTime) {
            console.log(patchTime + " " + patchBody.substring(1));
        }
        ip++;
        continue;
    }
    if (patchBody.startsWith("-")) {
        removeTime = patchTime;
        removeBodyList.push(patchBody.substring(1));
        ip++;
        continue;
    }
    break;
}



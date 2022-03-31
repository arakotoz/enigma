#!/usr/bin/env python3

import argparse
import os
import alienpy.alien as alien

""" 
---- How to run:
workon alice3.9.9
python3 getFileList.py -r 505713 -y 2021 -p NOV -u arakotoz
"""

def getUserHome(user):
    return "/alice/cern.ch/user/{}/{}".format(user[0], user)


def getDataDir(runNumber, period, year):
    return "/alice/data/{}/{}/{}/raw".format(year, period, runNumber)


def getRunFileList(runNumber, period, year, pattern, user, limit, isLocal):
    dataDir = getDataDir(runNumber, period, year)
    workDir = "{}/{}/{}".format(year, period, runNumber)
    outFilename = "{}".format(runNumber)
    if isLocal:
        outFilename += ".txt"
    else:
        outFilename += ".xml"
        workDir = "{}/mft/commissioning/".format(getUserHome(user)) + workDir

    return getFileList(dataDir, pattern, workDir, outFilename, limit, isLocal)


def getFileList(dataDir, pattern, workDir, outFilename, limit, isLocal):
    alien.setup_logging()
    jal = alien.AliEn()
    fullPath = "{}/{}".format(workDir, outFilename)
    findCmd = "find -l {}".format(limit)

    if isLocal:
        os.makedirs(workDir, exist_ok=True)
    else:
        findCmd += " -x {}".format(fullPath)
        out = jal.run("mkdir -p " + workDir)
        if (out.exitcode != 0):
            return out.exitcode
        out = jal.run("stat " + fullPath)
        if (out.exitcode == 0):
            # File exists
            jal.run("rm " + fullPath)
    findCmd += " {} {}".format(dataDir, pattern)
    print(findCmd)
    out = jal.run(findCmd)

    if isLocal:
        results = out.ansdict['results']
        print("Creating {} with {} entries".format(fullPath, len(results)))
        with open(fullPath, "w") as outFile:
            for dic in results:
                outFile.write("alien://{}\n".format(dic["lfn"]))

    return out.exitcode


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get file lists")
    parser.add_argument("--run", "-r", help="run number", required=True)
    parser.add_argument("--period", "-p", help="period", required=True)
    parser.add_argument("--year", "-y", help="year", default="2021")
    parser.add_argument("--pattern", "-f", help="user", default="o2*.root")
    parser.add_argument("--user", "-u", help="user", default="arakotoz")
    parser.add_argument(
        "--limit", "-l", help="maximum number of files", default="15000")
    parser.add_argument(
        "--local", help="produces a local text file", default=False, action="store_true")

    args = parser.parse_args()
    getRunFileList(args.run, args.period, args.year,
                   args.pattern, args.user, args.limit, args.local)

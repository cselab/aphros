#!/usr/bin/env python3

import argparse
import filecmp
import os
import shutil
import sys
import subprocess


class TestBase:
    """
    Base class for tests with reference data.
    """
    def __init__(self, description="", cases=[], verbose=True):
        """
        description: `str`
            Description shown in the help message.
        cases: `list(str)`
            Names of test cases.  If not empty, provides choices for `--run`.
        """

        self.verbose = verbose
        if self.verbose:

            def printlog(msg, end='\n'):
                sys.stderr.write(str(msg) + end)
                sys.stderr.flush()
        else:

            def printlog(msg, end='\n'):
                pass

        self.printlog = printlog

        description = " ".join(
            [description, "Default action is equivalent to '--run --check'."])
        parser = argparse.ArgumentParser(description=description)
        if cases:
            self.cases = cases
            parser.add_argument(
                '--run',
                action='store_true',
                help="Run selected test to generate output files."
                " Other options will use the test selected here.")
            parser.add_argument(
                'case',
                nargs='?',
                choices=cases,
                help="Command without any options provided implies"
                " '--run --check'")
        else:
            self.cases = None
            self.case = None
            parser.add_argument('--run',
                                action='store_true',
                                help="Run test to generate output files.")
        parser.add_argument('--check',
                            action='store_true',
                            help="Check output files against reference data."
                            " Exit with status 1 if failed.")
        parser.add_argument('--plot',
                            action='store_true',
                            help="Plot output data.")
        parser.add_argument('--plotref',
                            action='store_true',
                            help="Plot reference data.")
        parser.add_argument('--clean',
                            action='store_true',
                            help="Remove output files.")
        parser.add_argument('--update',
                            action='store_true',
                            help="Update reference data from output files.")
        parser.add_argument('--cases',
                            action='store_true',
                            help="Print names of all cases and exit.")
        parser.add_argument('--current',
                            action='store_true',
                            help="Print name of the current case and exit.")
        self.parser = parser

    def run(self, outdir, case=None):
        """
        Runs selected test case to generate output files.
        Returns a list of generated output files to include in reference data.
        """
        return []

    def check(self, outdir, refdir, output_files):
        """
        Checks output files against reference data.
        Returns False if failed.
        """
        r = True
        for f in output_files:
            out = os.path.join(outdir, f)
            ref = os.path.join(refdir, f)
            if not filecmp.cmp(out, ref):
                self.printlog("Files '{}' and '{}' differ".format(out, ref))
                r = False
        return r

    def update(self, outdir, refdir, output_files):
        """
        Updates reference data from output files.
        """
        os.makedirs(refdir, exist_ok=True)
        for f in output_files:
            self.copy_file(os.path.join(outdir, f), os.path.join(refdir, f))

    def plot(self, datadir, output_files):
        """
        Plots output files.
        datadir: `str`
            Directory with data files (output or reference).
        """
        self.printlog("plot: not implemented")
        pass

    def clean(self, outdir, output_files):
        """
        Removes output files.
        """
        ff = [os.path.join(outdir, f) for f in output_files] + ["testcase"]
        for f in ff:
            if os.path.isfile(f):
                os.remove(f)
                self.printlog("removed '{}'".format(f))

    def copy_file(self, src, dst):
        shutil.copy(src, dst)
        self.printlog("copied '{}' to '{}'".format(src, dst))

    def rename_file(self, src, dst):
        os.rename(src, dst)
        self.printlog("moved '{}' to '{}'".format(src, dst))

    def runcmd(self, cmd, echo=True):
        """
        Runs a command through the shell.
        cmd: `str`
            Command to run.
        """
        if echo:
            self.printlog(cmd)
        subprocess.run(cmd, shell=True, check=True)

    def _get_casedir(self, case, base):
        """
        Returns relative path to directory with case data.
        case: `str`
            Case name.
        base: `str`
            Base directory (e.g "output" or "ref").
        """
        if case:
            return os.path.join(base, case)
        return base

    def _move_to_outdir(self, outdir, output_files):
        """
        Moves output files to output directory.
        base: `str`
            Base directory (e.g "output" or "ref").
        """
        os.makedirs(outdir, exist_ok=True)
        for f in output_files:
            out = os.path.join(outdir, f)
            self.rename_file(f, out)

    def main(self):
        """
        Entry point. Parses command line arguments.
        If option --run, creates file 'testcase' with lines:
            - test case name
            - relative path to directory with output data
            - relative path to directory with reference data
            - names of files to include in reference data
        Example of 'testcase' file:
            vof
            ref/vof
            sm_0000.vtk
            traj_0000.vtk
        """
        self.args = self.parser.parse_args()
        args = self.args

        status = 0

        if args.cases:
            if self.cases:
                print(' '.join(self.cases))
            exit(status)

        casefile = "testcase"
        if args.current:
            if os.path.isfile(casefile):
                with open(casefile, 'r') as f:
                    print(f.readlines()[0])
            exit(status)

        if not any([
                args.run, args.check, args.plot, args.plotref, args.update,
                args.clean
        ]):
            if self.cases and args.case is None:
                self.parser.error(
                    "Running without options is equivalent to '--run --check',"
                    " but requres CASE as the only argument."
                    " Options are: {{{}}}".format(','.join(self.cases)))
            args.run = True
            args.check = True

        if args.run:
            if self.cases:
                self.case = args.case
                if self.case is None:
                    self.parser.error(
                        "Missing case to run. Options are {{{}}}".format(
                            ','.join(self.cases)))
                self.output_files = self.run(self.case)
            else:
                self.case = None
                self.output_files = self.run()
            self.outdir = self._get_casedir(self.case, "output")
            self.refdir = self._get_casedir(self.case, "ref")
            self._move_to_outdir(self.outdir, self.output_files)
            with open(casefile, 'w') as f:
                f.write(str(self.case) + '\n')
                f.write(self.outdir + '\n')
                f.write(self.refdir + '\n')
                for line in self.output_files:
                    f.write(line + '\n')
        else:
            assert os.path.isfile(casefile), \
                "File '{}' not found, use --run option first.".format(casefile)
            with open(casefile, 'r') as f:
                lines = f.readlines()
                lines = [l.strip() for l in lines]
                self.case = lines[0]
                self.outdir = lines[1]
                self.refdir = lines[2]
                self.output_files = lines[3:]
                if getattr(args, "case", None) is not None \
                        and args.case != self.case:
                    self.parser.error(
                        "Case '{}' differs from previously stored '{}'."
                        " Provide option --run to run the new case.".format(
                            args.case, self.case))

        if args.update:
            self.printlog("updating reference data in '{}'".format(
                self.refdir))
            self.update(self.outdir, self.refdir, self.output_files)

        if args.clean:
            self.clean(self.outdir, self.output_files)

        if args.plot:
            self.printlog("plotting output data in '{}'".format(self.outdir))
            self.plot(self.outdir, self.output_files)

        if args.plotref:
            self.printlog("plotting reference data in '{}'".format(
                self.refdir))
            self.plot(self.refdir, self.output_files)

        if args.check:
            if self.check(self.outdir, self.refdir, self.output_files):
                self.printlog("check with '{}' passed".format(self.refdir))
            else:
                self.printlog("check with '{}' failed".format(self.refdir))
                status = 1

        exit(status)

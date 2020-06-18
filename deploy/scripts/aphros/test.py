#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess


class TestBase:
    """
    Base class for tests with reference data.
    """
    def __init__(self, description="", cases=[]):
        """
        description: `str`
            Description shown in the help message.
        cases: `list(str)`
            Names of test cases.  If not empty, provides choices for `--run`.
        """
        description = " ".join(
            [description, "Default action is equivalent to '--run --check'."])
        parser = argparse.ArgumentParser(description=description)
        if cases:
            self.cases = cases
            parser.add_argument(
                '--run',
                action='store_true',
                help="Run selected test to generate output files."
                " Other commands will use the test selected here.")
            parser.add_argument(
                'case',
                nargs='?',
                choices=cases,
                help="If not other options provided, equivalent to"
                " '--run CASE --check'")
        else:
            self.cases = None
            parser.add_argument('--run',
                                action='store_true',
                                help="Run test to generate output files.")
        parser.add_argument('--check',
                            action='store_true',
                            help="Check output files against reference data."
                            " Exit with status 1 if failed.")
        parser.add_argument('--plot',
                            action='store_true',
                            help="Plot output files")
        parser.add_argument('--clean',
                            action='store_true',
                            help="Remove output files")
        parser.add_argument('--update',
                            action='store_true',
                            help="Update reference data from output files.")
        self.parser = parser

    def run(self, case=None):
        """
        Runs selected test case to generate output files.
        Returns a list of generated output files to include in reference data.
        """
        return []

    def check(self, refdir, output_files):
        """
        Checks output files against reference data.
        Returns False if failed.
        """
        r = True
        for out in output_files:
            ref = os.path.join(refdir, out)
            if not filecmp.cmp(out, ref):
                print("Files '{}' and '{}' differ".format(out, ref))
                r = False
        return r

    def update(self, refdir, output_files):
        """
        Updates reference data from output files.
        """
        os.makedirs(refdir, exist_ok=True)
        for out in output_files:
            ref = os.path.join(refdir, out)
            shutil.copy(out, ref)
            print("'{}' -> '{}'".format(out, ref))

    def plot(self, output_files):
        """
        Plots output files.
        """
        print("plot: not implemented")
        pass

    def clean(self, output_files):
        """
        Removes output files.
        """
        ff = output_files + ["testcase"]
        for out in ff:
            if os.path.isfile(out):
                os.remove(out)
                print("removed '{}'".format(out))

    def runcmd(self, cmd, echo=True):
        """
        Runs a command through the shell.
        cmd: `str`
            Command to run.
        """
        if echo:
            print(cmd)
        subprocess.run(cmd, shell=True, check=True)

    def _get_refdir(self):
        """
        Returns relative path to directory with reference data.
        Requires `self.cases` and `self.case`.
        """
        base = "ref"
        if self.cases:
            assert self.case in self.cases
            return os.path.join(base, self.case)
        return base

    def main(self):
        """
        Entry point. Parses command line arguments.
        If option --run, creates file 'testcase' with lines:
            - test case name
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
        if not any([args.run, args.check, args.plot, args.update, args.clean]):
            if self.cases and args.case is None:
                self.parser.error(
                    "Running without options is equivalent to '--run --check',"
                    " but requres CASE as the only argument."
                    " Options are: {{{}}}".format(','.join(self.cases)))
            args.run = True
            args.check = True

        casefile = "testcase"
        if args.run:
            if self.cases:
                self.case = args.case
                self.output_files = self.run(self.case)
                if self.case is None:
                    self.parser.error(
                        "Missing case to run. Options are {{{}}}".format(
                            ','.join(self.cases)))
            else:
                self.case = None
                self.output_files = self.run()
            self.refdir = self._get_refdir()
            with open(casefile, 'w') as f:
                f.write(str(self.case) + '\n')
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
                self.refdir = lines[1]
                self.output_files = lines[2:]
                if args.case is not None and args.case != self.case:
                    self.parser.error(
                        "Case '{}' differs from previously stored '{}'."
                        " Provide option --run to run the new case.".format(
                            args.case, self.case))

        if args.update:
            self.update(self.refdir, self.output_files)

        if args.clean:
            self.clean(self.output_files)

        if args.plot:
            self.plot(self.output_files)

        if args.check:
            if not self.check(self.refdir, self.output_files):
                status = 1

        exit(status)

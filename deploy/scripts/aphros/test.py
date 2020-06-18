#!/usr/bin/env python3

import argparse
import os
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
                choices=cases,
                help="Run selected test to generate output files."
                " Other commands will use the test selected here.")
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

    def check(self, refdir):
        """
        Checks output files against reference data.
        Returns False if failed.
        """
        return True

    def plot(self):
        """
        Plots output files.
        """
        pass

    def update(self):
        """
        Updates reference data from output files.
        """
        pass

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
        if not any([args.run, args.check, args.plot, args.update]):
            args.run = True
            args.check = True

        casefile = "testcase"
        if args.run:
            if len(self.cases):
                self.case = args.run
                self.output_files = self.run(self.case)
            else:
                self.case = None
                self.output_files = self.run()
            self.refdir = self._get_refdir()
            with open(casefile, 'w') as f:
                f.write(self.case + '\n')
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

        if args.check:
            if not self.check(self.refdir):
                status = 1

        if args.update:
            self.update()

        if args.plot:
            self.plot()

        exit(status)

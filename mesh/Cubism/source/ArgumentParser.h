/*
 *  ArgumentParser.h
 *  Cubism
 *
 *	This argument parser assumes that all arguments are optional ie, each of the argument names is preceded by a '-'
 *		all arguments are however NOT optional to avoid a mess with default values and returned values when not found!
 *
 *	More converter could be required:
 *		add as needed
 *			TypeName as{TypeName}() in Value
 *
 *  Created by Christian Conti on 6/7/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <limits>


using namespace std;

class Value
{
private:
	string content;

public:

	Value() : content("") {}

	Value(string content_) : content(content_) { /*printf("%s\n",content.c_str());*/ }

	double asDouble(double def=0)
	{
		if (content == "")
        {
            ostringstream sbuf;
            sbuf << def;
            content = sbuf.str();
        }
		return (double) atof(content.c_str());
	}

	int asInt(int def=0)
	{
		if (content == "")
        {
            ostringstream sbuf;
            sbuf << def;
            content = sbuf.str();
        }
		return atoi(content.c_str());
	}

	bool asBool(bool def=false)
	{
		if (content == "")
        {
            if (def) content = "true";
            else     content = "false";
        }
		if (content == "0") return false;
		if (content == "false") return false;

		return true;
	}

	string asString(string def="")
	{
		if (content == "") content = def;

		return content;
	}
};

class CommandlineParser
{
private:
	const int iArgC;
	const char** vArgV;
	bool bStrictMode, bVerbose;

protected:
	map<string,Value> mapArguments;

public:

	Value& operator()(const string arg)
	{
		if (bStrictMode)
		{
			map<string,Value>::const_iterator it = mapArguments.find(arg);

			if (it == mapArguments.end())
			{
				printf("Runtime option NOT SPECIFIED! ABORTING! name: %s\n",arg.data());
				abort();
			}
		}

		if (bVerbose) printf("%s is %s\n", arg.data(), mapArguments[arg].asString().data());
		return mapArguments[arg];
	}

	bool check(const string arg) const
	{
		return mapArguments.find(arg) != mapArguments.end();
	}

	CommandlineParser(const int argc, const char ** argv) : mapArguments(), iArgC(argc), vArgV(argv), bStrictMode(false), bVerbose(true)
	{
		for (int i=1; i<argc; i++)
			if (argv[i][0] == '-')
			{
				string values = "";
				int itemCount = 0;

				for (int j=i+1; j<argc; j++)
					if (argv[j][0] == '-')
						break;
					else
					{
						if (strcmp(values.c_str(), ""))
							values += ' ';

						values += argv[j];
						itemCount++;
					}

				if (itemCount == 0)
					values = "true";

				mapArguments[argv[i]] = Value(values);
				i += itemCount;
			}

		mute();
		//printf("found %ld arguments of %d\n",mapArguments.size(),argc);
	}

	int getargc() const { return iArgC; }

	const char** getargv() const { return vArgV; }

	void set_strict_mode()
	{
		bStrictMode = true;
	}

	void unset_strict_mode()
	{
		bStrictMode = false;
	}

	void mute()
	{
		bVerbose = false;
	}

	void loud()
	{
		bVerbose = true;
	}

	void save_options(string path=".")
	{
		string options;
		for(map<string,Value>::iterator it=mapArguments.begin(); it!=mapArguments.end(); it++)
		{
			options+= it->first + " " + it->second.asString() + " ";
		}
		string filepath = (path + "/" + string("argumentparser.log"));
		FILE * f = fopen(filepath.data(), "a");
		if (f == NULL)
		{
			printf("impossible to write %s.\n", filepath.data());
			return;
		}
		fprintf(f, "%s\n", options.data());
		fclose(f);
	}

	void print_args()
	{
		for(map<string,Value>::iterator it=mapArguments.begin(); it!=mapArguments.end(); it++)
		{
            std::cout.width(50);
            std::cout.fill('.');
            std::cout << std::left << it->first;
            std::cout << ": " << it->second.asString() << std::endl;
		}
    }
};


class ArgumentParser: public CommandlineParser
{
    typedef std::map<std::string, Value> ArgMap;
    typedef std::map<std::string, Value*> pArgMap;
    typedef std::map<std::string, ArgMap* > FileMap;

    // keep a reference form option origin
    ArgMap  from_commandline;
    FileMap from_files;
    pArgMap from_code;

    // helper
    void _ignoreComments(std::ifstream& stream, const char commentChar = '#')
    {
        stream >> std::ws;
        int nextchar = stream.peek();
        while (nextchar == commentChar)
        {
            stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            stream >> std::ws;
            nextchar = stream.peek();
        }
    }

    bool _existKey(std::string& key) const
    {
        const std::string og_key = key;
        bool bExist = true;

        // look for both possible keys (i.e. with leading "-" and without)
        ArgMap::const_iterator it = mapArguments.find(key);
        if (it == mapArguments.end())
        {
            if (key[0] == '-') key = key.erase(0, 1);
            else key = "-" + key;
            it = mapArguments.find(key);
            if (it == mapArguments.end())
            {
                key = og_key;
                bExist = false;
            }
        }
        return bExist;
    }

    inline std::string _stripKey(std::string key) const
    {
        if (key[0] == '-') key = key.erase(0, 1);
        return key;
    }


public:
    ArgumentParser(const int _argc, const char ** _argv):
        CommandlineParser(_argc, _argv) { from_commandline = mapArguments; }

    virtual ~ArgumentParser()
    {
        for (FileMap::iterator it = from_files.begin(); it != from_files.end(); it++)
            delete it->second;
    }

    void readFile(const std::string filepath)
    {
        from_files[filepath] = new ArgMap;
        ArgMap& myFMap = *(from_files[filepath]);

        std::ifstream confFile(filepath.c_str());
        if (confFile.is_open())
        {
            // read (key value) pairs from input file, ignore comments
            // beginning with "#"
            std::string key, val;
            while (!confFile.eof())
            {
                _ignoreComments(confFile);
                confFile >> key >> val;
                if (_existKey(key)) continue;
                std::pair<string, Value> item(key, Value(val));
                mapArguments.insert(item); // add to parent container
                myFMap.insert(item); // add to private container
            }
        }
    }

    Value& operator()(std::string key)
    {
        const bool bDefaultInCode = !_existKey(key);
        Value& retval = CommandlineParser::operator()(key);
        if (bDefaultInCode) from_code[key] = &retval;
        return retval;
    }

    inline bool check(std::string key) const { return _existKey(key); }
    inline bool exist(std::string key) const { return _existKey(key); }

    void print_args()
    {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "* Summary:" << std::endl;
        std::cout << "*    Parameter read from command line:                " << from_commandline.size() << std::endl;
        size_t nFiles = 0;
        size_t nFileParameter = 0;
        for (FileMap::const_iterator it=from_files.begin(); it!=from_files.end(); ++it)
        {
            if (it->second->size() > 0)
            {
                ++nFiles;
                nFileParameter += it->second->size();
            }
        }
        std::cout << "*    Parameter read from " << std::setw(3) << std::right << nFiles << " file(s):                 " << nFileParameter << std::endl;
        std::cout << "*    Parameter read from defaults in code:            " << from_code.size() << std::endl;
        std::cout << "*    Total number of parameter read from all sources: " << mapArguments.size() << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

        // command line given arguments
        if (!from_commandline.empty())
        {
            std::cout << "* Command Line:" << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            for(ArgMap::iterator it=from_commandline.begin(); it!=from_commandline.end(); it++)
            {
                std::cout.width(50);
                std::cout.fill('.');
                std::cout << std::left << _stripKey(it->first);
                std::cout << ": " << it->second.asString() << std::endl;
            }
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        }

        // options read from input files
        if (!from_files.empty())
        {
            for (FileMap::iterator itFile=from_files.begin(); itFile!=from_files.end(); itFile++)
            {
                if (!itFile->second->empty())
                {
                    std::cout << "* File: " << itFile->first << std::endl;
                    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                    ArgMap& fileArgs = *(itFile->second);
                    for(ArgMap::iterator it=fileArgs.begin(); it!=fileArgs.end(); it++)
                    {
                        std::cout.width(50);
                        std::cout.fill('.');
                        std::cout << std::left << _stripKey(it->first);
                        std::cout << ": " << it->second.asString() << std::endl;
                    }
                    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                }
            }
        }

        // defaults defined in code
        if (!from_code.empty())
        {
            std::cout << "* Defined in Code:" << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            for(pArgMap::iterator it=from_code.begin(); it!=from_code.end(); it++)
            {
                std::cout.width(50);
                std::cout.fill('.');
                std::cout << std::left << _stripKey(it->first);
                std::cout << ": " << it->second->asString() << std::endl;
            }
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        }
    }
};

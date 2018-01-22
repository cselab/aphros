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
#include <ctime>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <limits>


using namespace std;

class Value
{
private:
	string content;

public:

	Value() : content("") {}

	Value(string content_) : content(content_) { /*printf("%s\n",content.c_str());*/ }

    Value(const Value& c) : content(c.content) {}

    Value& operator=(const Value& rhs)
    {
        if (this != &rhs)
            content = rhs.content;
        return *this;
    }
    Value& operator+=(const Value& rhs)
    {
        content += " " + rhs.content;
        return *this;
    }
    Value operator+(const Value& rhs) { return Value(content + " " + rhs.content); }

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

    friend std::ostream& operator<<(std::ostream& lhs, const Value& rhs)
    {
        lhs << rhs.content;
        return lhs;
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

    inline void _normalizeKey(std::string& key) const
    {
        if (key[0] == '-') key.erase(0,1);
        if (key[0] == '+') key.erase(0,1);
    }

    inline bool _existKey(const std::string& key, const std::map<std::string,Value>& container) const
    {
        std::map<std::string,Value>::const_iterator it = container.find(key);
        return it != container.end();
    }

public:

	Value& operator()(string key)
	{
        _normalizeKey(key);
		if (bStrictMode)
		{
			if (!_existKey(key,mapArguments))
			{
				printf("Runtime option NOT SPECIFIED! ABORTING! name: %s\n",key.data());
				abort();
			}
		}

		if (bVerbose) printf("%s is %s\n", key.data(), mapArguments[key].asString().data());
		return mapArguments[key];
	}

	inline bool check(string key) const
	{
        _normalizeKey(key);
		return _existKey(key,mapArguments);
	}

	CommandlineParser(const int argc, const char ** argv) : mapArguments(), iArgC(argc), vArgV(argv), bStrictMode(false), bVerbose(true)
	{
		for (int i=1; i<argc; i++)
			if (argv[i][0] == '-')
			{
				string values = "";
				int itemCount = 0;

				for (int j=i+1; j<argc; j++)
                {
                    const bool leadingDash = (argv[j][0] == '-');
                    const char c = argv[j][1];
                    const bool firstNumeric = ((c >= '0' && c <= '9') || c == 0) ? true : false;
					if (leadingDash && !firstNumeric)
						break;
					else
					{
						if (strcmp(values.c_str(), ""))
							values += ' ';

						values += argv[j];
						itemCount++;
					}
                }

				if (itemCount == 0)
					values = "true";

                std::string key(argv[i]);
                key.erase(0,1); // remove leading '-'
                if (key[0] == '+')
                {
                    key.erase(0,1);
                    if (!_existKey(key,mapArguments))
                        mapArguments[key] = Value(values); // skip leading white space
                    else
                        mapArguments[key] += Value(values);
                }
                else
                {
                    if (!_existKey(key,mapArguments))
                        mapArguments[key] = Value(values);
                }

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

    const char commentStart;

    // keep a reference from option origin
    ArgMap  from_commandline;
    FileMap from_files;
    pArgMap from_code;

    // for runtime interaction (we keep the original map)
    ArgMap mapRuntime;

    // helper
    void _ignoreComments(std::istream& stream, const char commentChar)
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

    void _parseFile(std::ifstream& stream, ArgMap& container)
    {
        // read (key value) pairs from input file, ignore comments
        // beginning with commentStart
        _ignoreComments(stream, commentStart);
        while (!stream.eof())
        {
            std::string line, key, val;
            std::getline(stream, line);
            std::istringstream lineStream(line);
            lineStream >> key;
            lineStream >> val;
            _ignoreComments(lineStream, commentStart);
            while(!lineStream.eof())
            {
                std::string multiVal;
                lineStream >> multiVal;
                val += (" " + multiVal);
                _ignoreComments(lineStream, commentStart);
            }

            const Value V(val);
            if (key[0] == '-')
                key.erase(0,1);

            if (key[0] == '+')
            {
                key.erase(0,1);
                if (!_existKey(key,container)) // skip leading white space
                    container[key] = V;
                else
                    container[key] += V;
            }
            else if (!_existKey(key,container))
                container[key] = V;
            _ignoreComments(stream, commentStart);
        }
    }


public:
    ArgumentParser(const int _argc, const char ** _argv, const char cstart='#'):
        CommandlineParser(_argc, _argv), commentStart(cstart)
    {
        from_commandline = mapArguments;
    }

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
        if (confFile.good())
        {
            _parseFile(confFile, mapArguments);
            confFile.clear();
            confFile.seekg(0, ios::beg);
            _parseFile(confFile, myFMap); // we keep a reference for each separate file read
        }
        confFile.close();
    }

    Value& operator()(std::string key)
    {
        _normalizeKey(key);
        const bool bDefaultInCode = !_existKey(key,mapArguments);
        Value& retval = CommandlineParser::operator()(key);
        if (bDefaultInCode) from_code[key] = &retval;
        return retval;
    }

    inline bool exist(std::string key) const { return check(key); }

    void write_runtime_environment() const
    {
        time_t rawtime;
        std::time(&rawtime);
        struct tm* timeinfo = std::localtime(&rawtime);
        char buf[256];
        std::strftime(buf, 256, "%A, %h %d %Y, %r", timeinfo);

        std::ofstream runtime("runtime_environment.conf");
        runtime << commentStart << " RUNTIME ENVIRONMENT SETTINGS" << std::endl;
        runtime << commentStart << " ============================" << std::endl;
        runtime << commentStart << " " << buf << std::endl;
        runtime << commentStart << " Use this file to set runtime parameter interactively." << std::endl;
        runtime << commentStart << " The parameter are read every \"refreshperiod\" steps." << std::endl;
        runtime << commentStart << " When editing this file, you may use comments and string concatenation." << std::endl;
        runtime << commentStart << " The simulation can be terminated without killing it by setting \"exit\" to true." << std::endl;
        runtime << commentStart << " (This will write a serialized restart state. Set \"exitsave\" to false if not desired.)" << std::endl;
        runtime << commentStart << std::endl;
        runtime << commentStart << " !!! WARNING !!! EDITING THIS FILE CAN POTENTIALLY CRASH YOUR SIMULATION !!! WARNING !!!" << std::endl;
        for (typename std::map<std::string,Value>::const_iterator it = mapArguments.begin(); it != mapArguments.end(); ++it)
            runtime << it->first << '\t' << it->second << std::endl;
    }

    void read_runtime_environment()
    {
        mapRuntime.clear();
        std::ifstream runtime("runtime_environment.conf");
        if (runtime.good())
            _parseFile(runtime, mapRuntime);
        runtime.close();
    }

    Value& parseRuntime(std::string key)
    {
        _normalizeKey(key);
        if (!_existKey(key,mapRuntime))
        {
            printf("ERROR: Runtime parsing for key %s NOT FOUND!! Check your runtime_environment.conf file\n",key.data());
            abort();
        }
        return mapRuntime[key];
    }

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
                std::cout << std::left << it->first;
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
                        std::cout << std::left << it->first;
                        std::cout << ": " << it->second.asString() << std::endl;
                    }
                    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                }
            }
        }

        // defaults defined in code
        if (!from_code.empty())
        {
            std::cout << "* Defaults in Code:" << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            for(pArgMap::iterator it=from_code.begin(); it!=from_code.end(); it++)
            {
                std::cout.width(50);
                std::cout.fill('.');
                std::cout << std::left << it->first;
                std::cout << ": " << it->second->asString() << std::endl;
            }
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        }
    }
};

/*
 *  Profiler.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <assert.h>
#undef min
#undef max
#include <vector>
#undef min
#undef max
#include <map>
#include <string>
#include <stdio.h>
#include <stack>

using namespace std;

#include <sys/time.h>
//#include <tbb/tick_count.h>
//namespace tbb { class tick_count; }

const bool bVerboseProfiling = false;

class ProfileAgent
{
//	typedef tbb::tick_count ClockTime;
	typedef timeval ClockTime;
	
	enum ProfileAgentState{ ProfileAgentState_Created, ProfileAgentState_Started, ProfileAgentState_Stopped};

	ClockTime m_tStart, m_tEnd;
	ProfileAgentState m_state;
	long double m_dAccumulatedTime;
	int m_nMeasurements;
	int m_nMoney;

	static void _getTime(ClockTime& time);

	static double _getElapsedTime(const ClockTime& tS, const ClockTime& tE);

	void _reset()
	{
		m_tStart = ClockTime();
		m_tEnd = ClockTime();
		m_dAccumulatedTime = 0;
		m_nMeasurements = 0;
		m_nMoney = 0;
		m_state = ProfileAgentState_Created;
	}

public:

	ProfileAgent():m_tStart(), m_tEnd(), m_state(ProfileAgentState_Created),
	m_dAccumulatedTime(0), m_nMeasurements(0), m_nMoney(0) {}

	void start()
	{
		assert(m_state == ProfileAgentState_Created || m_state == ProfileAgentState_Stopped);

		if (bVerboseProfiling) {printf("start\n");}

		_getTime(m_tStart);

		m_state = ProfileAgentState_Started;
	}

	void stop(int nMoney=0)
	{
		assert(m_state == ProfileAgentState_Started);

		if (bVerboseProfiling) {printf("stop\n");}

		_getTime(m_tEnd);
		m_dAccumulatedTime += _getElapsedTime(m_tStart, m_tEnd);
		m_nMeasurements++;
		m_nMoney += nMoney;
		m_state = ProfileAgentState_Stopped;
	}

	friend class Profiler;
};

struct ProfileSummaryItem
{
	string sName;
	double dTime;
	double dAverageTime;
	int nMoney;
	int nSamples;

	ProfileSummaryItem(string sName_, double dTime_, int nMoney_, int nSamples_):
		sName(sName_), dTime(dTime_), nMoney(nMoney_),nSamples(nSamples_), dAverageTime(dTime_/nSamples_){}
};


class Profiler
{
protected:

	map<string, ProfileAgent*> m_mapAgents;
	stack<string> m_mapStoppedAgents;
		
public:
	void push_start(string sAgentName)
	{
		if (m_mapStoppedAgents.size() > 0)
			getAgent(m_mapStoppedAgents.top()).stop();
			
		m_mapStoppedAgents.push(sAgentName);
		getAgent(sAgentName).start();
	}
	
	void pop_stop()
	{
		string sCurrentAgentName = m_mapStoppedAgents.top();
		getAgent(sCurrentAgentName).stop();
		m_mapStoppedAgents.pop();
		
		if (m_mapStoppedAgents.size() == 0) return;
		
		getAgent(m_mapStoppedAgents.top()).start();
	}
	
	void clear()
	{
		for(map<string, ProfileAgent*>::iterator it = m_mapAgents.begin(); it != m_mapAgents.end(); it++)
		{
			delete it->second;

			it->second = NULL;
		}

		m_mapAgents.clear();
	}

	Profiler(): m_mapAgents(){}

	~Profiler()
	{
		clear();
	}

	void printSummary(FILE *outFile=NULL) const
	{
		vector<ProfileSummaryItem> v = createSummary();

		double dTotalTime = 0;
		double dTotalTime2 = 0;
		for(vector<ProfileSummaryItem>::const_iterator it = v.begin(); it!= v.end(); it++)
			dTotalTime += it->dTime;

		for(vector<ProfileSummaryItem>::const_iterator it = v.begin(); it!= v.end(); it++)
		dTotalTime2 += it->dTime - it->nSamples*1.30e-6;

		for(vector<ProfileSummaryItem>::const_iterator it = v.begin(); it!= v.end(); it++)
		{
			const ProfileSummaryItem& item = *it;
			const double avgTime = item.dAverageTime;

			printf("[%15s]: \t%02.0f-%02.0f%%\t%03.3e (%03.3e) s\t%03.3f (%03.3f) s\t(%d samples)\n",
				   item.sName.data(), 100*item.dTime/dTotalTime, 100*(item.dTime- item.nSamples*1.3e-6)/dTotalTime2, avgTime,avgTime-1.30e-6,  item.dTime, item.dTime- item.nSamples*1.30e-6, item.nSamples);
			if (outFile) fprintf(outFile,"[%15s]: \t%02.2f%%\t%03.3f s\t(%d samples)\n",

				   item.sName.data(), 100*item.dTime/dTotalTime, avgTime, item.nSamples);
		}

		printf("[Total time]: \t%f\n", dTotalTime);
		if (outFile) fprintf(outFile,"[Total time]: \t%f\n", dTotalTime);
		if (outFile) fflush(outFile);
		if (outFile) fclose(outFile);
	}

	vector<ProfileSummaryItem> createSummary(bool bSkipIrrelevantEntries=true) const
	{
		vector<ProfileSummaryItem> result;
		result.reserve(m_mapAgents.size());

		for(map<string, ProfileAgent*>::const_iterator it = m_mapAgents.begin(); it != m_mapAgents.end(); it++)
		{
			const ProfileAgent& agent = *it->second;
			if (!bSkipIrrelevantEntries || agent.m_dAccumulatedTime>1e-3)
				result.push_back(ProfileSummaryItem(it->first, agent.m_dAccumulatedTime, agent.m_nMoney, agent.m_nMeasurements));
		}

		return result;
	}

	void reset()
	{
		printf("reset\n");
		for(map<string, ProfileAgent*>::const_iterator it = m_mapAgents.begin(); it != m_mapAgents.end(); it++)
			it->second->_reset();
	}

	ProfileAgent& getAgent(string sName)
	{
		if (bVerboseProfiling) {printf("%s ", sName.data());}

		map<string, ProfileAgent*>::const_iterator it = m_mapAgents.find(sName);

		const bool bFound = it != m_mapAgents.end();

		if (bFound) return *it->second;

		ProfileAgent * agent = new ProfileAgent();

		m_mapAgents[sName] = agent;

		return *agent;
	}

	friend class ProfileAgent;
};




/*
 *  OutputProcessing.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 09/12/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef OUTPUTPROCESSING_H_CRQLO4OA
#define OUTPUTPROCESSING_H_CRQLO4OA

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Tests.h"
#include "ArgumentParser.h"
#include "Types.h"
#include "Profiler.h"
#include "Streamer.h"
#include "BlockProcessor_OMP.h"

// dumper (we do not support wavelet compressed output w/o MPI)
#include "HDF5Dumper.h"
#include "HDF5SliceDumper.h"


#define __REGISTER_H5_ENTITY__(KEY, INFO, TYPE, HEAVY, TENTITY, TDUMP, TGRID, ENTITY, GRID, PARS, FDUMP, PROCESSINGELEMENT, PROC, LAB) \
    this->_register(Item(#KEY"_a2",      new ProcessingElement<TENTITY,TDUMP,TGRID>(#INFO": Alpha 2", &ENTITY, &GRID, &PARS, &FDUMP<PROCESSINGELEMENT, StreamerAlpha2>, TYPE, HEAVY)));


template <typename TGrid>
struct Slice
{
    typedef TGrid GridType;

    TGrid * grid;
    int id;
    int dir;
    int idx;
    int width, height;
    bool valid;
    Slice() : grid(NULL), id(-1), dir(-1), idx(-1), width(0), height(0), valid(false) {}

    template <typename TSlice>
    static std::vector<TSlice> getSlices(ArgumentParser& parser, TGrid& grid)
    {
        typedef typename TGrid::BlockType B;
        int Dim[3];
        Dim[0] = grid.getBlocksPerDimension(0)*B::sizeX;
        Dim[1] = grid.getBlocksPerDimension(1)*B::sizeY;
        Dim[2] = grid.getBlocksPerDimension(2)*B::sizeZ;

        std::vector<TSlice> slices(0);
        const size_t nSlices = parser("nslices").asInt(0);
        for (size_t i = 0; i < nSlices; ++i)
        {
            TSlice thisOne;
            thisOne.id = i+1;
            thisOne.grid = &grid;
            assert(thisOne.grid != NULL);

            std::ostringstream identifier;
            identifier << "slice" << i+1;
            // fetch direction
            const std::string sDir = identifier.str() + "_direction";
            if (parser.check(sDir)) thisOne.dir = parser(sDir).asInt(0);
            const bool bDirOK = (thisOne.dir >= 0 && thisOne.dir < 3);
            assert(bDirOK);

            // compute index
            const std::string sIndex = identifier.str() + "_index";
            const std::string sFrac  = identifier.str() + "_fraction";
            if (parser.check(sIndex)) thisOne.idx = parser(sIndex).asInt(0);
            else if (parser.check(sFrac))
            {
                const double fraction = parser(sFrac).asDouble(0.0);
                const int idx = static_cast<int>(Dim[thisOne.dir] * fraction);
                thisOne.idx = (fraction == 1.0) ? Dim[thisOne.dir]-1 : idx;
            }
            const bool bIdxOK = (thisOne.idx >= 0 && thisOne.idx < Dim[thisOne.dir]);
            assert(bIdxOK);

            if (bDirOK && bIdxOK) thisOne.valid = true;
            else
            {
                std::cerr << "Slice: WARNING: Ill defined slice \"" << identifier.str() << "\"... Skipping this one" << std::endl;
                thisOne.valid = false;
                slices.push_back(thisOne);
                continue;
            }

            // define slice layout
            if (thisOne.dir == 0)
            {
                thisOne.width  = Dim[2];
                thisOne.height = Dim[1];
            }
            else if (thisOne.dir == 1)
            {
                thisOne.width  = Dim[2];
                thisOne.height = Dim[0];
            }
            else if (thisOne.dir == 2)
            {
                thisOne.width  = Dim[0];
                thisOne.height = Dim[1];
            }
            slices.push_back(thisOne);
        }
        return slices;
    }
};


template <typename TSlice>
class SliceProcessor
{
    typedef void (*Tkern)(const TSlice&, const int, const Real, const std::string, const std::string, const bool);

public:
    SliceProcessor(ArgumentParser& parser, typename TSlice::GridType& grid, const bool verbose=true, Tkern f = NULL) :
        m_verbose(verbose), m_writer(f)
    {
        m_slices = TSlice::template getSlices<TSlice>(parser, grid);
    }
    ~SliceProcessor() {}

    inline void setFunctor(Tkern f)
    {
        m_writer = f;
    }

    inline void process(const int sliceID, const int stepID, const Real t, const std::string fname, const std::string dpath=".", const bool bXMF=true)
    {
        // TODO: (fabianw@mavt.ethz.ch; Thu 29 Sep 2016 02:59:54 PM CEST)
        // WARNING: sliceID starts from zero, eventhough slice numbering in the
        // .conf file starts from 1 (to be consistent with earlier
        // implementation of sensors)
        if (sliceID >= 0 && sliceID < static_cast<int>(m_slices.size()))
            _process(m_slices[sliceID], stepID, t, fname, dpath, bXMF);
        else
            if (m_verbose)
                std::cerr << "SliceProcessor: WARNING: Requesting undefined slice... Skipping slice" << sliceID << std::endl;
    }

    inline void process_all(const int stepID, const Real t, const std::string fname, const std::string dpath=".", const bool bXMF=true)
    {
        for (size_t i = 0; i < m_slices.size(); ++i)
            _process(m_slices[i], stepID, t, fname, dpath, bXMF);
    }

    void showSlices()
    {
        if (m_verbose)
        {
            std::cout << "Got n = " << m_slices.size() << " slices:" << std::endl;
            for (size_t i = 0; i < m_slices.size(); ++i)
            {
                std::cout << "Slice ID = " << m_slices[i].id << std::endl;
                std::cout << "\tSlice index:     " << m_slices[i].idx << std::endl;
                std::cout << "\tSlice direction: " << m_slices[i].dir << std::endl;
                std::cout << "\tSlice dimension: (" << m_slices[i].width << ", " << m_slices[i].height << ")" << std::endl;
                std::cout << "\tSlice valid:     " << m_slices[i].valid << std::endl;
                std::cout << std::endl;
            }
        }
    }

private:
    const bool m_verbose;
    Tkern m_writer;
    std::vector<TSlice> m_slices;

    inline void _process(const TSlice& slice, const int stepID, const Real t, const std::string fname, const std::string path, const bool bXMF)
    {
        if (m_writer && slice.valid)
            m_writer(slice, stepID, t, fname, path, bXMF);
        if (m_verbose && !m_writer)
            std::cerr << "SliceProcessor: WARNING: No functor defined... Skipping slice" << slice.id << std::endl;
    }
};


template <typename TGrid>
class OperatorTypeBase
{
public:
    OperatorTypeBase(const std::string name, const int c, const unsigned int mask, const unsigned int inval) :
        m_name(name), m_class(c), m_mask(mask), m_invalid(inval)
    { }
    virtual ~OperatorTypeBase() {}

    inline const std::string& name() const { return m_name; }
    inline int operatorClass() const { return m_class; }
    inline unsigned int mask() const { return m_mask; }
    inline unsigned int invalidate() const { return m_invalid; }
    virtual void operator()(TGrid& grid, const Real t, const bool verbose) = 0;

private:
    const std::string m_name;
    const int m_class;
    const unsigned int m_mask;
    const unsigned int m_invalid;
};


template <typename TGrid, typename TKernel>
class OperatorType : public OperatorTypeBase<TGrid>
{
    typedef void (*TFunctor)(TKernel, TGrid&, const Real, const bool);

public:
    OperatorType(TFunctor f) :
        OperatorTypeBase<TGrid>(TKernel::NAME,TKernel::CLASS,TKernel::MASK,TKernel::INVALID),
        m_functor(f)
    { }
    virtual ~OperatorType() {}

    virtual void operator()(TGrid& grid, const Real t, const bool verbose)
    {
        m_functor(m_kernel, grid, t, verbose);
    }

private:
    TFunctor m_functor;
    TKernel m_kernel;
};


class ProcessingElementBase
{
public:
    enum EType { H5=0, H5SLICE, VP };

    ProcessingElementBase(const std::string& name, const EType type, const bool heavy=false, const int operatorClass=0) :
    m_name(name), m_operatorClass(operatorClass), m_heavy(heavy), m_type(type), m_dispatch(true)
    { }
    virtual ~ProcessingElementBase() {}

    inline const std::string& name() const { return m_name; }
    inline int operatorClass() const { return m_operatorClass; }
    inline bool is_heavy() const { return m_heavy; }
    inline void set_dispatch(const bool v=true) { m_dispatch = v; }
    inline bool dispatch() const { return m_dispatch; }
    inline EType type() const { return m_type; }

    // interface
    virtual void operator()(const Real t, const bool verbose) = 0;
    virtual void operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool bXMF) = 0;
    virtual void dispose() = 0;
    virtual std::string operatorName() const = 0;
    virtual unsigned int operatorInvalidate() const = 0;
    virtual unsigned int operatorMask() const = 0;

private:
    const std::string m_name;
    const int m_operatorClass;
    const bool m_heavy;
    const EType m_type;

protected:
    bool m_dispatch;
};


template <typename TEntity, typename TOutFunctor, typename TGrid>
class ProcessingElement : public ProcessingElementBase
{
public:
    ProcessingElement(const std::string& name, TEntity* e, TGrid* g, ArgumentParser* p, TOutFunctor f,
            const EType type, const bool heavy, const int oClass=0, OperatorTypeBase<TGrid>* op=NULL) :
    ProcessingElementBase(name, type, heavy, oClass), m_entity(e), m_outputFunctor(f), m_grid(g), m_operator(op), m_parser(p)
    { }
    virtual ~ProcessingElement() { if (m_operator) delete m_operator; }

    virtual void operator()(const Real t, const bool verbose)
    {
        if (m_operator)
        {
            assert(this->operatorClass() == m_operator->operatorClass());
            m_operator->operator()(*m_grid, t, verbose);
        }
    }

    virtual void operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
    {
        std::cout << "ERROR: ProcessingElement: operator()(...) needs specialization!" << std::endl;
        abort();
    }

    virtual void dispose() { }
    virtual std::string operatorName() const
    {
        if (m_operator)
            return m_operator->name();
        return std::string("Requires no operator");
    }
    virtual unsigned int operatorInvalidate() const
    {
        if (m_operator)
            return m_operator->invalidate();
        return ~0x0;
    }
    virtual unsigned int operatorMask() const
    {
        if (m_operator)
            return m_operator->mask();
        return 0x0;
    }

private:

    TEntity* const m_entity;
    TOutFunctor m_outputFunctor;
    TGrid* const m_grid;
    OperatorTypeBase<TGrid>* m_operator;
    ArgumentParser* const m_parser;
};

typedef Slice<Grid_t> SliceGrid_t;
typedef void (*Dumper_h5_full_grid)(const Grid_t&, const int, const Real, const std::string, const std::string, const bool);
typedef void (*Dumper_h5_slice_grid)(const SliceGrid_t&, const int, const Real, const std::string, const std::string, const bool);

// specialization
template <>
void ProcessingElement<Grid_t, Dumper_h5_full_grid, Grid_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_outputFunctor(*m_grid, step_id, t, basename, path, true);
}

template <>
void ProcessingElement<SliceProcessor<SliceGrid_t>, Dumper_h5_slice_grid, Grid_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_entity->setFunctor(m_outputFunctor);
    m_entity->process_all(step_id, t, basename, path);
}


typedef ProcessingElementBase PElement;
typedef std::pair<std::string, PElement*> Item;
typedef std::map<std::string, PElement*> TMap;
typedef TMap::iterator TMapItem; // Item * const
typedef TMap::const_iterator TMapConstItem; // const Item * const

bool pipe_sort_order(TMapConstItem a, TMapConstItem b)
{
    const int ca = a->second->operatorClass();
    const int cb = b->second->operatorClass();
    if (ca == cb)
        return a->second->operatorName() < b->second->operatorName();
    return ca < cb;
}

template <typename TGrid, template <typename> class TSlice=Slice>
class OutputProcessing
{
public:
    OutputProcessing(ArgumentParser& p, TGrid& grid, const bool verbose=true) :
        m_registered(false), m_parser(p), m_sliceProcessor(p,grid,verbose),
        m_bVerbose(verbose), m_bBypass(false),
        m_dumpcount(0), m_steplast(0), m_timelast(0.0)
    { }

    virtual ~OutputProcessing()
    {
        if (m_registered)
        {
            for (typename TMap::iterator it = m_items.begin(); it != m_items.end(); ++it)
            {
                it->second->dispose();
                delete it->second;
            }
        }

        if (m_bVerbose)
            m_log.close();
    }

    inline void toggle_bypass() { m_bBypass = 1 - m_bBypass; }
    inline void list_registered() const
    {
        if (m_bVerbose) _list_all();
    }
    inline void register_all(TGrid& grid)
    {
        if (!m_registered)
        {
            _set_parameter();
            _register_all(grid);

            if (m_bVerbose)
            {
                m_log.open("output_processor.log", std::ofstream::app);
                m_log << "#                                 DUMP LOG                                     " << std::endl;
                m_log << "# =============================================================================" << std::endl;
            }

            m_registered = true;
        }
    }

    Real operator()(const int step_id, const Real t, const int maxsteps, const Real tend, Profiler& prof, const bool bInvalidateOP=true, const bool bVeto=false)
    {
        const Real dt_step = fabs(tend - t);
        const Real dt_next = fabs(m_tdump - t);
        Real dt = (dt_next < dt_step) ? dt_next : dt_step;
        if (m_bIO)
        {
            const bool bExtrema    = (step_id == 0 || step_id == maxsteps);
            const bool bDumpByStep = bExtrema || (m_dumpperiod != 0 && step_id % m_dumpperiod == 0);
            const bool bDumpByTime = (dt_next <= std::numeric_limits<Real>::epsilon()) || (dt_step <= std::numeric_limits<Real>::epsilon());
            if ((bDumpByStep || bDumpByTime) && !bVeto || m_bBypass)
            {
                // concept:
                // (1) prepare processing pipe (sorted by operator class)
                // (2) some infos
                // (3) compute all output element wise.  m_validOperators contains
                // all valid operators evaluated in _process_all() which can be
                // re-used in later parts of this step_id (if needed).
                // (4) update m_steplast m_timelast

                // 1.)
                const bool bFinal         = (step_id == maxsteps) || (dt_step <= std::numeric_limits<Real>::epsilon());
                const bool bDispatchHeavy = (m_heavySkipStep > 0 ? (m_dumpcount % m_heavySkipStep == 0) : true);
                const bool bDispatchAll   = bDispatchHeavy || bFinal || m_bBypass;
                _prepare(bDispatchAll);

                // 2.)
                if (m_bVerbose)
                {
                    std::cout << "GENERATING OUTPUT:" << std::endl;
                    _to_stream(std::cout, step_id, t, bDumpByStep, bDumpByTime, bDispatchAll);

                    m_log << "Dump count: " << std::setfill('0') << std::setw(6) << std::right << m_dumpcount << std::endl;
                    _to_stream(m_log, step_id, t, bDumpByStep, bDumpByTime, bDispatchAll);
                    m_log << std::endl;
                }

                // 3.)
                const std::string path = m_parser("fpath").asString();
                std::ostringstream basename;
                basename << "data_" << std::setfill('0') << std::setw(6) << std::right << step_id;
                prof.push_start("I/O");
                _process_all(step_id, t, basename.str(), path);
                prof.pop_stop();

                // 4.)
                m_steplast = step_id;
                m_timelast = t;
                ++m_dumpcount;

                if (m_bVerbose)
                    std::cout << "GENERATING OUTPUT: DONE" << std::endl;
            }

            if (bDumpByTime)
            {
                m_tdump = std::min(m_tdump + m_dumpdt, tend);
                const Real dt_next = fabs(m_tdump - t);
                dt = (dt_next < dt_step) ? dt_next : dt_step;
            }
        }
        return dt;
    }


private:

    template <typename, typename, template <typename> class>
    friend class Test_SteadyState;

    bool m_registered;

    inline void _list_all() const
    {
        for (TMapConstItem it = m_items.begin(); it != m_items.end(); ++it)
            std::cout << it->second->name() << " (key = " << it->first << ")" << std::endl;
    }


protected:

    // class local types
    typedef TSlice<TGrid> USlice;
    typedef void (*TFunc_h5_grid)(const TGrid&, const int, const Real, const std::string, const std::string, const bool);
    typedef void (*TFunc_h5_slice)(const USlice&, const int, const Real, const std::string, const std::string, const bool);
    typedef SliceProcessor<USlice> TSP;

    TMap m_items;
    ArgumentParser& m_parser;
    TSP m_sliceProcessor;
    const bool m_bVerbose;
    bool m_bBypass; // bypass all flags (force full dump)
    std::ofstream m_log;
    std::vector<TMapItem> m_processingPipe;

    // these attributes are also modified through friends
    std::string m_channels;
    int m_dumpperiod, m_dumpcount, m_steplast;
    int m_heavySkipStep; // skip resource and memory heavy dumps every n steps
    Real m_tdump, m_dumpdt, m_timelast;
    bool m_bIO, m_bVP, m_bHDF, m_bHDF_SLICE;
    std::pair<int, std::vector<std::string> > m_validOperators;

    // helper
    inline void _register(Item item) { m_items.insert(item); }

    virtual void _register_h5_grid(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5, HDF5_Full, PElement::H5, true, TGrid, TFunc_h5_grid, TGrid, grid, grid, m_parser, DumpHDF5, TGrid, processOMP, Lab)
    }

    virtual void _register_h5_slice(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5s, HDF5_Slice, PElement::H5SLICE, false, TSP, TFunc_h5_slice, TGrid, m_sliceProcessor, grid, m_parser, DumpSliceHDF5, USlice, processOMP, Lab)
    }

    virtual void _register_all(TGrid& grid)
    {
        _register_h5_grid(grid);
        _register_h5_slice(grid);
    }

    void _process_all(const int step_id, const Real t, const std::string& basename, const std::string& path)
    {
        for (size_t i=0; i<m_processingPipe.size(); ++i)
        {
            PElement* const item = m_processingPipe[i]->second;
            if (item->dispatch())
            {
                // concept:
                // (1) check if operator evaluation is required
                // (2) process output

                // 1.)
                if (item->operatorClass() > 0)
                {
                    const int c = item->operatorClass();
                    int& invalidate = m_validOperators.first;
                    invalidate &= item->operatorInvalidate();
                    invalidate |= item->operatorMask();

                    std::vector<std::string>& lastOfClass = m_validOperators.second;
                    if (item->operatorName() != lastOfClass[c])
                    {
#ifndef NDEBUG
                        if (m_bVerbose)
                            std::cout << "* Processing operator: " << item->operatorName() << std::endl;
#endif
                        (*item)(t, m_bVerbose);
                        lastOfClass[c] = item->operatorName();
                    }
                }

                // 2.)
                if (m_bVerbose)
                    std::cout << "* Processing: " << item->name() << std::endl;
                (*item)(step_id, t, basename, path, m_bVerbose);
            }
        }
    }

    void _prepare(const bool bDispatchAll)
    {
        // concept:
        // (1): split channels into vector
        // (2): assign items into processing pipe
        // (3): sort pipe based on operator class (evaluate in ascending class
        // order)
        // (4): prepare items for dispatch

        // 1.)
        std::vector<std::string> definedChannels(0);
        std::istringstream iss(m_channels);
        iss >> std::ws;
        while (!iss.eof())
        {
            std::string channel;
            iss >> channel;
            definedChannels.push_back(channel);
            iss >> std::ws;
        }

        // 2.)
        m_processingPipe.clear();
        m_processingPipe.reserve(definedChannels.size());
        for (typename std::vector<std::string>::const_iterator channel = definedChannels.begin(); channel != definedChannels.end(); ++channel)
        {
            TMapItem item = m_items.find(*channel);
            if (item == m_items.end()) continue;
            m_processingPipe.push_back(item);
        }

        // 3.)
        std::sort(m_processingPipe.begin(), m_processingPipe.end(), pipe_sort_order);

        // 4.)
        for (size_t i=0; i<m_processingPipe.size(); ++i)
        {
            PElement * const item = m_processingPipe[i]->second;
            const bool dispatch = !(item->is_heavy()) || bDispatchAll;
            item->set_dispatch(dispatch);
            if (item->type() == PElement::H5)
                item->set_dispatch(dispatch && m_bHDF);
            else if (item->type() == PElement::H5SLICE)
                item->set_dispatch(dispatch && m_bHDF_SLICE);
            else if (item->type() == PElement::VP)
                item->set_dispatch(dispatch && m_bVP);
        }
    }

    inline void _set_parameter()
    {
        m_dumpperiod    = m_parser("dumpperiod").asInt(0);
        m_dumpdt        = m_parser("dumpdt").asDouble(m_parser("tend").asDouble());
        m_tdump         = m_parser("tdump").asDouble(m_dumpdt);
        m_heavySkipStep = m_parser("heavyskipstep").asInt(1);
        m_bIO           = m_parser("io").asBool(false);
        m_bVP           = m_parser("vp").asBool(false);
        m_bHDF          = m_parser("hdf").asBool(false);
        m_bHDF_SLICE    = m_parser("hdf_slice").asBool(false);
        m_channels      = m_parser("channels").asString("none");
    }

    inline void _to_stream(std::ostream& stream, const int step_id, const Real t, const bool bDumpByStep, const bool bDumpByTime, const bool bDispatchAll)
    {
        const std::string dispatch[2] = {"hold", "dispatch"};

        stream << '\t' << m_processingPipe.size() << " Items";
        if (0 == m_processingPipe.size())
            stream << std::endl;
        else
        {
            stream << " =" << std::endl;
            for (size_t i = 0; i < m_processingPipe.size(); ++i)
                stream << "\t\t" << std::setfill('.') << std::setw(24) << std::left << m_processingPipe[i]->first << "(state = " << dispatch[m_processingPipe[i]->second->dispatch()] << ")" << std::endl;
        }
        stream << "\tStep = " << step_id << " (" << step_id-m_steplast << " steps passed since last dump)" << std::endl;
        stream << "\tTime = " << t << " (" << t-m_timelast << " time passed since last dump)" << std::endl;
        stream << "\tbDumpByStep  = " << bDumpByStep << std::endl;
        stream << "\tbDumpByTime  = " << bDumpByTime << std::endl;
        stream << "\tbDispatchAll = " << bDispatchAll << std::endl;
        stream << "\tbBypass      = " << m_bBypass << std::endl;
    }
};

#endif /* OUTPUTPROCESSING_H_CRQLO4OA */

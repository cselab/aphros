/*
 *  Test_Statistics.h
 *  CUBISMTests
 *
 *  Created by Jonas Sukys on May 10, 2015.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Types.h"
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#include "Reader_WaveletCompression.h"
#include <ArgumentParser.h>
#include <GridMPI.h>
#define VERBOSE 1

typedef GridMPI< FluidGrid > G;


class Test_Statistics : public Simulation
{
protected:

	int BPDX, BPDY, BPDZ;
	int XPESIZE, YPESIZE, ZPESIZE;
  double GAMMA1, GAMMA2;
	int step_id;
  double t;
  double dt;
	bool VERBOSITY;
	G * grid;
	ArgumentParser parser;
	string inputfile_name, outputfile_name;
  string inputpath, outputpath;
	int myrank;
  
  Real extent;
  Real extents [3];
  
  bool BC_PERIODIC [3];
  
  // TODO: add all required streamers
	SerializerIO_WaveletCompression_MPI_SimpleBlocking <G, StreamerGridPointIterative> mywaveletdumper;
  
	void _ic(FluidGrid& grid)
	{
    const int NCHANNELS = 7;
    
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		
		MPI_Comm comm  = MPI_COMM_WORLD;
		MPI::Intracomm& mycomm = MPI::COMM_WORLD;
		const bool swapbytes = parser.check("-swap");
    
    // read all channels
    Reader_WaveletCompressionMPI* vpreader [NCHANNELS];
    for (int channel = 0; channel < NCHANNELS; channel++) {
      
      std::stringstream file;
      file << inputpath << "/datawavelet";
      file.setf (ios::dec | ios::right);
      file.width(6);
      file.fill('0');
      file << step_id;
      
      std::stringstream suffix;
      suffix << ".GridPointIterative.channel" << channel;
      
      string fullname = file.str() + suffix.str();
		  
      cout << "inputfile_name: "  << fullname << "\n" ;
      
  		vpreader [channel] = new Reader_WaveletCompressionMPI (mycomm, fullname, swapbytes);
      vpreader [channel] -> load_file();
    }
    
		static Real targetdata [NCHANNELS] [_BLOCKSIZE_] [_BLOCKSIZE_] [_BLOCKSIZE_];
    
		//#pragma omp parallel for
		for (int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
      
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      
      for (int channel = 0; channel < NCHANNELS; channel++)
        vpreader [channel] -> load_block2 (info.index[0], info.index[1], info.index[2], targetdata [channel]);
			
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
					{
						b(ix, iy, iz).r = targetdata [0] [iz][iy][ix];
            b(ix, iy, iz).u = targetdata [1] [iz][iy][ix];
            b(ix, iy, iz).v = targetdata [2] [iz][iy][ix];
            b(ix, iy, iz).w = targetdata [3] [iz][iy][ix];
            b(ix, iy, iz).p = targetdata [4] [iz][iy][ix];
            b(ix, iy, iz).G = targetdata [5] [iz][iy][ix];
            b(ix, iy, iz).P = targetdata [6] [iz][iy][ix];
					}
		}
    
		// cleanup
    for (int channel = 0; channel < NCHANNELS; channel++)
      delete vpreader [channel];
	}


	void _setup_mpi_constants(int& xpesize, int& ypesize, int& zpesize)
	{
		xpesize = parser("-xpesize").asInt(1);
		ypesize = parser("-ypesize").asInt(1);
		zpesize = parser("-zpesize").asInt(1);
    
    const int BPD_PE_MAX = max (max (BPDX * xpesize, BPDY * ypesize), BPDZ * zpesize);
    extents [0] = extent * (BPDX * xpesize) / BPD_PE_MAX;
    extents [1] = extent * (BPDY * ypesize) / BPD_PE_MAX;
    extents [2] = extent * (BPDZ * zpesize) / BPD_PE_MAX;
	}
    
public:
	const bool isroot;
	
	Test_Statistics (const bool isroot, const int argc, const char ** argv) : isroot(isroot) , grid(NULL), parser(argc,argv) { }
    

	void setup()
	{
		parser.unset_strict_mode();

		BPDX = parser("-bpdx").asInt(1);
		BPDY = parser("-bpdy").asInt(1);
		BPDZ = parser("-bpdz").asInt(1);
    
    GAMMA1 = parser("-g1").asDouble(6.59);
    GAMMA2 = parser("-g2").asDouble(1.4);
    
    step_id = parser("-stepid").asInt(-1);
    
    t = parser("-t").asDouble(0.0);
    
    dt = 0;
    
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
		inputpath = parser("-inputpath").asString("none");
		outputpath = parser("-outputpath").asString("none");
    
    extent = parser("-extent").asDouble(1.0);
    
    const int BPD_MAX = max (max (BPDX, BPDY), BPDZ);
    extents [0] = extent * BPDX / (double) BPD_MAX;
    extents [1] = extent * BPDY / (double) BPD_MAX;
    extents [2] = extent * BPDZ / (double) BPD_MAX;
    
    BC_PERIODIC [0] = parser("-bcperiodicx").asInt(0);
    BC_PERIODIC [1] = parser("-bcperiodicy").asInt(0);
    BC_PERIODIC [2] = parser("-bcperiodicz").asInt(0);
    
		if ((inputfile_name == "none")||(outputfile_name == "none") || step_id == -1)
		{
			printf("usage: %s -inputpath <path> -outputpath <path> -stepid <number> [-swap]\n", "tool");
			exit(1);
		}
    
		_setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);
    
		if (!isroot)
			VERBOSITY = 0;
    
		if (VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("////////////       TEST Statistics MPI       ///////////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}

		grid = new G (XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);

		assert(grid != NULL);
    
    //const int NCHANNELS = 7;
    //for (int channel = 0; channel < NCHANNELS, channel++)
    //  _ic (*grid, channel);
    _ic (*grid);
    
    dumpStatistics (*grid, step_id, t, dt);
    
		vp (*grid, step_id);
	}
  
  // structure for the "quantity of interest"
  // consecutive fields 'value' and 'rank' ensure proper MPI communication
  struct QoI {
    string name;
    double value;
    int    rank;
    double pos[3];
    double distance;
    double distance_x;
    double distance_y;
    double distance_z;
    MPI::Op op;
    
    // constructor for integral quantities
    QoI (const string &name, const double &value) : name(name), value(value), op(MPI::SUM) {}
    
    // constructor for extremal quantities
    QoI (const string &name, const double &value, const MPI::Op &op) : name(name), value(value), op(op) {
      rank = MPI::COMM_WORLD.Get_rank();
    }
    
    // update minimum
    inline void update_min (const double &new_value, const Real (&new_pos) [3]) {
      if (new_value < value) {
        value = new_value;
        pos[0] = (double)new_pos[0];
        pos[1] = (double)new_pos[1];
        pos[2] = (double)new_pos[2];
      }
    }
    
    // update maximum
    inline void update_max (const double &new_value, const Real (&new_pos) [3]) {
      if (new_value > value) {
        value = new_value;
        pos[0] = (double)new_pos[0];
        pos[1] = (double)new_pos[1];
        pos[2] = (double)new_pos[2];
      }
    }
  };
  
  void dumpStatistics (G& grid, const int step_id, const Real t, const Real dt)
  {
    // vector of pointers to integral quantities
    vector <QoI*> qi;
    
    // integral quantities (position is not relevant)
    QoI r_avg  ("r_avg",  0.); qi.push_back (&r_avg);  // density
    QoI r2_avg ("r2_avg", 0.); qi.push_back (&r2_avg); // density of phase 2
    QoI u_avg  ("u_avg",  0.); qi.push_back (&u_avg);
    QoI v_avg  ("v_avg",  0.); qi.push_back (&v_avg);
    QoI w_avg  ("w_avg",  0.); qi.push_back (&w_avg);
    QoI m_avg  ("m_avg",  0.); qi.push_back (&m_avg);  // magnitude of velocity vector
    //QoI W_avg  ("W_avg",  0.); qi.push_back (&W_avg);  // magnitude of vorticity vector
    QoI ke_avg ("ke_avg", 0.); qi.push_back (&ke_avg);
    QoI e_avg  ("e_avg",  0.); qi.push_back (&e_avg);
    QoI p_avg  ("p_avg",  0.); qi.push_back (&p_avg);
    QoI pw_avg ("pw_avg", 0.); qi.push_back (&pw_avg);
    QoI c_avg  ("c_avg",  0.); qi.push_back (&c_avg);
    QoI M_avg  ("M_avg",  0.); qi.push_back (&M_avg);
    QoI V2     ("V2",     0.); qi.push_back (&V2);     // volume of phase 2
    QoI Req    ("Req",    0.); qi.push_back (&Req);    // equivalent radius of phase 2
    QoI Vc     ("Vc",     0.); qi.push_back (&Vc);     // volume of the whole cloud (phases 1 and 2)
    
    // vector of pointers to extremal quantities
    vector <QoI*> qe;
    
    // extremal quantities (position is also relevant)
    QoI r_min  ("r_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&r_min);
    QoI r_max  ("r_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&r_max);
    QoI r2_min ("r2_min",  HUGE_VAL, MPI::MINLOC); qe.push_back (&r2_min);
    QoI r2_max ("r2_max", -HUGE_VAL, MPI::MAXLOC); qe.push_back (&r2_max);
    QoI u_min  ("u_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&u_min);
    QoI u_max  ("u_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&u_max);
    QoI v_min  ("v_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&v_min);
    QoI v_max  ("v_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&v_max);
    QoI w_min  ("w_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&w_min);
    QoI w_max  ("w_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&w_max);
    QoI m_min  ("m_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&m_min);
    QoI m_max  ("m_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&m_max);
    QoI W_min  ("W_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&W_min);
    //QoI W_max  ("W_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&W_max);
    QoI ke_min ("ke_min",  HUGE_VAL, MPI::MINLOC); qe.push_back (&ke_min);
    QoI ke_max ("ke_max", -HUGE_VAL, MPI::MAXLOC); qe.push_back (&ke_max);
    QoI e_min  ("e_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&e_min);
    QoI e_max  ("e_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&e_max);
    QoI p_min  ("p_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&p_min);
    QoI p_max  ("p_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&p_max);
    QoI c_min  ("c_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&c_min);
    QoI c_max  ("c_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&c_max);
    QoI M_min  ("M_min",   HUGE_VAL, MPI::MINLOC); qe.push_back (&M_min);
    QoI M_max  ("M_max",  -HUGE_VAL, MPI::MAXLOC); qe.push_back (&M_max);
    QoI pw_min ("pw_min",  HUGE_VAL, MPI::MINLOC); qe.push_back (&pw_min);
    QoI pw_max ("pw_max", -HUGE_VAL, MPI::MAXLOC); qe.push_back (&pw_max);
    
    vector<BlockInfo> g_vInfo = grid.getBlocksInfo();
    const double h = g_vInfo.front().h_gridpoint;
    const double h3 = h*h*h;
    const double V_block = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ;
    
    double V2_block;
    Real pos [3];
    double r, u, v, w, p, G, P;
    double W0, W1, W2;
    double r2, m, W, ke, e, c, M;
    int phase;
    
    int blocks = 0;
    int wall_cells = 0;
    
    for(int i=0; i<(int)g_vInfo.size(); i++)
    {
      BlockInfo info = g_vInfo[i];
      
      bool bBoundary = 0;
      for (int d = 0; d < 3; d++)
        if ((info.index[d]==0 && ! BC_PERIODIC[d]) || (info.index[d]==grid.getBlocksPerDimension(d)-1 && ! BC_PERIODIC[d]))
          bBoundary = 1;
      
      parser.unset_strict_mode();
      // ignore sponge blocks
      if(parser("-sponge").asInt(0)==1 && bBoundary)
        continue;
      
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      blocks++;
      V2_block = 0;
      
      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
        for(int iy=0; iy<FluidBlock::sizeY; iy++)
          for(int ix=0; ix<FluidBlock::sizeX; ix++)
          {
            r = (double)b(ix, iy, iz).r;
            u = (double)b(ix, iy, iz).u;
            v = (double)b(ix, iy, iz).v;
            w = (double)b(ix, iy, iz).w;
            p = (double)b(ix, iy, iz).p;
            G = (double)b(ix, iy, iz).G;
            P = (double)b(ix, iy, iz).P;
            
            phase = G > (0.1*min(1/(GAMMA1-1),1/(GAMMA2-1))+0.9*max(1/(GAMMA1-1),1/(GAMMA2-1)))? 2 : 1;
            
            /*
            // vorticity components
            W0 = (double)b.tmp[ix][iy][iz][0];
            W1 = (double)b.tmp[ix][iy][iz][1];
            W2 = (double)b.tmp[ix][iy][iz][2];
            */
            
            // === update integral quantities
            
            r_avg.value += r;
            u_avg.value += u;
            v_avg.value += v;
            w_avg.value += w;
            p_avg.value += p;
            
            // TODO: this is unclear
            r2 = r * (1-min(max((G-1/(GAMMA1-2))/(1/(GAMMA1-1)-1/(GAMMA2-1)),(double)0),(double)1));
            r2_avg.value += r2;
            
            m = sqrt (u*u + v*v + w*w);
            m_avg.value += m;
            
            /*
            W = sqrt (W0*W0 + W1*W1 + W2*W2);
            W_avg.value += W;
            */
            
            ke = 0.5 * r * (u*u + v*v + w*w);
            ke_avg.value += ke;
            
            e = ke + G * p + P;
            e_avg.value += e;
            
            c = sqrt((1/G + 1)*(p + P/G/(1/G + 1))/r);
            c_avg.value += c;
            
            M = m / c;
            M_avg.value += M;
            
            if (phase == 2) V2_block += 1;
            
            // === update extremal quantities
            
            // get position
            info.pos (pos, ix, iy, iz);
            
            r_min.update_min (r, pos);
            r_max.update_max (r, pos);
            r2_min.update_min (r2, pos);
            r2_max.update_max (r2, pos);
            u_min.update_min (u, pos);
            u_max.update_max (u, pos);
            v_min.update_min (v, pos);
            v_max.update_max (v, pos);
            w_min.update_min (w, pos);
            w_max.update_max (w, pos);
            ke_min.update_min (ke, pos);
            ke_max.update_max (ke, pos);
            m_min.update_min (m, pos);
            m_max.update_max (m, pos);
            //W_min.update_min (W, pos);
            //W_max.update_max (W, pos);
            e_min.update_min (e, pos);
            e_max.update_max (e, pos);
            p_min.update_min (p, pos);
            p_max.update_max (p, pos);
            c_min.update_min (c, pos);
            c_max.update_max (c, pos);
            M_min.update_min (M, pos);
            M_max.update_max (M, pos);
            
            if (info.index[2]==0 && iz==0) {
              pw_avg.value += p;
              wall_cells++;
              pw_min.update_min (p, pos);
              pw_max.update_max (p, pos);
            }
          }
      
      V2.value += V2_block;
      if (V2_block > 1e-3 * V_block)
        Vc.value += V_block;
    }
    
    int rank = MPI::COMM_WORLD.Get_rank();
    int root = rank == 0;
    
    // reduce integral quantities
    for (int i = 0; i < qi.size(); i++)
      MPI::COMM_WORLD.Reduce (root ? MPI::IN_PLACE : &(qi[i]->value), &(qi[i]->value), 1, MPI::DOUBLE, qi[i]->op, 0);
    
    /*
     // MPI datatype for packing value and rank
     const int    nitems = 2;
     int          blocklengths [2] = {1, 1};
     MPI::Datatype types [2] = {MPI::DOUBLE, MPI::INT};
     MPI::Datatype mpi_value_rank;
     MPI::Aint     offsets [2];
     offsets [0] = 0;
     offsets [1] = (char*)&(qe[0]->rank) - (char*)&(qe[0]->value);
     mpi_value_rank = MPI::Datatype::Create_struct (nitems, blocklengths, offsets, types);
     mpi_value_rank.Commit ();
     */
    
    // reduce extremal quantities
    for (int i = 0; i < qe.size(); i++)
      MPI::COMM_WORLD.Allreduce (MPI::IN_PLACE, &(qe[i]->value), 1, MPI::DOUBLE_INT, qe[i]->op);
    
    // destruct mpi_QoI
    //mpi_value_rank.Free ();
    
    // get positions of extremal quantities
    MPI::Status status;
    for (int i = 0; i < qe.size(); i++)
      if ( qe[i]->rank != 0 ) {
        if ( rank == qe[i]->rank )
          MPI::COMM_WORLD.Send (&(qe[i]->pos), 3, MPI::DOUBLE, 0, i);
        else if ( root )
          MPI::COMM_WORLD.Recv (&(qe[i]->pos), 3, MPI::DOUBLE, qe[i]->rank, i, status);
      }
    
    // compute distances from the center of the domain
    double center [3];
    center[0] = 0.5 * extents[0];
    center[1] = 0.5 * extents[1];
    center[2] = 0.5 * extents[2];
    for (int i = 0; i < qe.size(); i++) {
      qe[i]->distance = sqrt ( (center[0] - qe[i]->pos[0]) * (center[0] - qe[i]->pos[0]) + (center[1] - qe[i]->pos[1]) * (center[1] - qe[i]->pos[1]) + (center[2] - qe[i]->pos[2]) * (center[2] - qe[i]->pos[2]) );
      qe[i]->distance_x = sqrt ( (center[1] - qe[i]->pos[1]) * (center[1] - qe[i]->pos[1]) + (center[2] - qe[i]->pos[2]) * (center[2] - qe[i]->pos[2]) );
      qe[i]->distance_y = sqrt ( (center[0] - qe[i]->pos[0]) * (center[0] - qe[i]->pos[0]) + (center[2] - qe[i]->pos[2]) * (center[2] - qe[i]->pos[2]) );
      qe[i]->distance_z = sqrt ( (center[0] - qe[i]->pos[0]) * (center[0] - qe[i]->pos[0]) + (center[1] - qe[i]->pos[1]) * (center[1] - qe[i]->pos[1]) );
    }
    
    // reduce the total number of blocks
    MPI::COMM_WORLD.Reduce (root ? MPI::IN_PLACE : &blocks, &blocks, 1, MPI::DOUBLE, MPI::SUM, 0);
    
    // reduce the total number of wall cells
    MPI::COMM_WORLD.Reduce (root ? MPI::IN_PLACE : &wall_cells, &wall_cells, 1, MPI::DOUBLE, MPI::SUM, 0);
    
    // finalize intergral quantities
    double factor = 1.0 / (blocks * FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ);
    r_avg.value  *= factor;
    r2_avg.value *= factor;
    u_avg.value  *= factor;
    v_avg.value  *= factor;
    w_avg.value  *= factor;
    m_avg.value  *= factor;
    //W_avg.value  *= factor;
    ke_avg.value *= factor;
    e_avg.value  *= factor;
    p_avg.value  *= factor;
    if (wall_cells > 0)
      pw_avg.value /= wall_cells;
    c_avg.value  *= factor;
    M_avg.value  *= factor;
    V2.value     *= h3;
    Vc.value     *= h3;
    Req.value     = pow (0.75 * V2.value / M_PI, 1./3.);
    
    // write output file
    if (root)
    {
      string filename = "statistics_pp.dat";
      
      // write header
      std::ofstream f;
      f.open (filename.c_str(), std::ofstream::in);
      if ( f.good() )
        f.close();
      else {
        f.close();
        f.open (filename.c_str(), std::ofstream::out);
        f << "step t dt";
        for (int i = 0; i < qi.size(); i++)
          f << " " << qi[i]->name;
        for (int i = 0; i < qe.size(); i++) {
          f << " " << qe[i]->name;
          f << " " << qe[i]->name << "_pos_x";
          f << " " << qe[i]->name << "_pos_y";
          f << " " << qe[i]->name << "_pos_z";
          f << " " << qe[i]->name << "_pos_d";
          f << " " << qe[i]->name << "_pos_d_x";
          f << " " << qe[i]->name << "_pos_d_y";
          f << " " << qe[i]->name << "_pos_d_z";
        }
        f << std::endl;
        f.close();
      }
      
      f.open (filename.c_str(), std::ofstream::app);
      f << step_id << " " << t << " " << dt;
      f << std::scientific;
      for (int i = 0; i < qi.size(); i++)
        f << " " << qi[i]->value;
      for (int i = 0; i < qe.size(); i++) {
        f << " " << qe[i]->value;
        f << " " << qe[i]->pos[0];
        f << " " << qe[i]->pos[1];
        f << " " << qe[i]->pos[2];
        f << " " << qe[i]->distance;
        f << " " << qe[i]->distance_x;
        f << " " << qe[i]->distance_y;
        f << " " << qe[i]->distance_z;
      }
      f << std::endl;
      f.close();
    }
  }
  
	void vp(G& grid, const int step_id)
	{
    if (isroot) cout << "dumping MPI VP is disabled for now\n" ;
    return;
    // TODO: copy vp routine from Cubism
		
    if (isroot) cout << "dumping MPI VP ...\n" ;

		const string path = parser("-fpath").asString(".");

		std::stringstream streamer;
		streamer<<path;
		streamer<<"/";
		streamer<<outputfile_name;
		streamer.setf(ios::dec | ios::right);
		streamer.width(6);
		streamer.fill('0');
		streamer<<step_id;

		double threshold = parser("-threshold").asDouble(1e-5);

		mywaveletdumper.verbose();
		printf("setting threshold to %f\n", threshold);
		mywaveletdumper.set_threshold(threshold);
		mywaveletdumper.Write<0>(grid, streamer.str());

		if (isroot) cout << "done" << endl;
	}

	void dispose()
	{
		if (grid!=NULL) {
			delete grid;
			grid = NULL;
		}
	}
};

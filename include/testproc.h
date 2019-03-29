// testproc.h

#ifndef testproc_h
#define testproc_h 1

#include "lcfiplus.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TH2D.h"

namespace lcfiplus {

class ZHHAlgo : public Algorithm {
 public:
  struct data {
    int mchdecaypdg[2];
    int mchbb;
    int mcnb;

    double thrust;
    double thaxis[3];
    double ycuts[10];
    int ntr;
    int npfo;

    double mass[15];
    double ntrjetmin;
    double pmiss[3];
    double emiss;

    double bcat[6];
    double btag[6];
    double ctag[6];
    double ejet[6];
    double pxjet[6];
    double pyjet[6];
    double pzjet[6];
    double ntrjet[6];
    double twovtxprobjet[6];
    double vtxangle[6];
    double mcnb6[6];
    double mcnc6[6];

    double bcat4[4];
    double btag4[4];
    double ctag4[4];
    double ejet4[4];
    double pxjet4[4];
    double pyjet4[4];
    double pzjet4[4];
    double ntrjet4[4];
    double twovtxprobjet4[4];
    double vtxangle4[4];

    double bcat5[5];
    double btag5[5];
    double ctag5[5];
    double ejet5[5];
    double pxjet5[5];
    double pyjet5[5];
    double pzjet5[5];
    double ntrjet5[5];
    double twovtxprobjet5[5];
    double vtxangle5[5];

    double bcat7[7];
    double btag7[7];
    double ctag7[7];
    double ejet7[7];
    double pxjet7[7];
    double pyjet7[7];
    double pzjet7[7];
    double ntrjet7[7];
    double twovtxprobjet7[7];
    double vtxangle7[7];

    double bcat8[8];
    double btag8[8];
    double ctag8[8];
    double ejet8[8];
    double pxjet8[8];
    double pyjet8[8];
    double pzjet8[8];
    double ntrjet8[8];
    double twovtxprobjet8[8];
    double vtxangle8[8];

    double bcatnv4[4];
    double btagnv4[4];
    double ctagnv4[4];
    double ejetnv4[4];
    double pxjetnv4[4];
    double pyjetnv4[4];
    double pzjetnv4[4];
    double ntrjetnv4[4];
    double twovtxprobjetnv4[4];
    double vtxanglenv4[4];

    double bcatnv5[5];
    double btagnv5[5];
    double ctagnv5[5];
    double ejetnv5[5];
    double pxjetnv5[5];
    double pyjetnv5[5];
    double pzjetnv5[5];
    double ntrjetnv5[5];
    double twovtxprobjetnv5[5];
    double vtxanglenv5[5];

    double bcatnv6[6];
    double btagnv6[6];
    double ctagnv6[6];
    double ejetnv6[6];
    double pxjetnv6[6];
    double pyjetnv6[6];
    double pzjetnv6[6];
    double ntrjetnv6[6];
    double twovtxprobjetnv6[6];
    double vtxanglenv6[6];

    double bcatnv7[7];
    double btagnv7[7];
    double ctagnv7[7];
    double ejetnv7[7];
    double pxjetnv7[7];
    double pyjetnv7[7];
    double pzjetnv7[7];
    double ntrjetnv7[7];
    double twovtxprobjetnv7[7];
    double vtxanglenv7[7];

    double bcatnv8[8];
    double btagnv8[8];
    double ctagnv8[8];
    double ejetnv8[8];
    double pxjetnv8[8];
    double pyjetnv8[8];
    double pzjetnv8[8];
    double ntrjetnv8[8];
    double twovtxprobjetnv8[8];
    double vtxanglenv8[8];

  };

  ZHHAlgo() {}
  virtual ~ZHHAlgo() {}

  static bool sortBtag(const Jet* j1, const Jet* j2) {
    double btag1 = j1->getParam("lcfiplus")->get<double>("BTag");
    double btag2 = j2->getParam("lcfiplus")->get<double>("BTag");

    return btag1 > btag2;
  }

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(ZHHAlgo,1);

 private:
  TFile* _file;

  string _jetname;
  string _jetname4;
  string _jetname5;
  string _jetname7;
  string _jetname8;

  string _jetnamenv4;
  string _jetnamenv5;
  string _jetnamenv6;
  string _jetnamenv7;
  string _jetnamenv8;

  JetVec* _jets;  //!
  JetVec* _jets4;  //!
  JetVec* _jets5;  //!
  JetVec* _jets7;  //!
  JetVec* _jets8;  //!

  JetVec* _jetsnv4;  //!
  JetVec* _jetsnv5;  //!
  JetVec* _jetsnv6;  //!
  JetVec* _jetsnv7;  //!
  JetVec* _jetsnv8;  //!

  TTree* _tree;

  data _d;

};

class TestAlgo : public Algorithm {
 public:
  TestAlgo() {}
  virtual ~TestAlgo() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(TestAlgo,1);

 private:
  TNtupleD* _nt;
  int _nev;

  TNtuple* _ntJet2;
  TNtuple* _nbJet;
  TFile* _file;

  string _v0vtxname;
  string _privtxname;
  string _jetname;
  bool _bbhh;

  VertexVec* _vertices;  //!
  VertexVec* _v0vertices;  //!
  JetVec* _jets;  //!

  // for old version

  string _vtxname;
  int _vtxsel;
  int _refine;
  TNtupleD* _ntJet;

  TH2D* _h;
  TH2D* _he;

};

class VertexAnalysis : public Algorithm {
 public:
  VertexAnalysis() {}
  virtual ~VertexAnalysis() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(VertexAnalysis,1);

 private:
  TNtupleD* _nt;
  int _nev;

  TFile* _file;

  string _privtxname;
  string _secvtxname;
};

class FlavtagReader : public Algorithm {
 public:
  FlavtagReader() {}
  virtual ~FlavtagReader() {}

  void init(Parameters* param);
  void process();
  void end();

  ClassDef(FlavtagReader,1);

 private:
  TNtupleD* _nt;
  TNtupleD* _ntev;
  int _nev;

  TNtuple* _ntJet2;
  TNtuple* _nbJet;
  TFile* _file;

  string _v0vtxname;
  string _privtxname;
  string _jetname;
  bool _bbhh;

  VertexVec* _vertices;  //!
  VertexVec* _v0vertices;  //!
  JetVec* _jets;  //!

  // for old version

  string _vtxname;
  int _vtxsel;
  int _refine;
  TNtupleD* _ntJet;

  TH2D* _h;
  TH2D* _he;

};

class TestAlgoV0 : public Algorithm {
 public:
  TestAlgoV0() {}
  virtual ~TestAlgoV0() {}
  void init(Parameters* param);
  void process();
  void end();
  ClassDef(TestAlgoV0,1);
 private:
  TTree* _ntp;
  TFile* _file;
  VertexVec* _vertices;  //!
  string _vtxname;

  struct VtxData {
    double x;
    double y;
    double z;
    double r;
    double cs;
    double phi;
    double chrg;
    double dirdot;
    double dirdot2;
    int ntrk;
    double mks;
    double ml0;
    double mconv;
    double mks2;
    double ml02;
    int v0;
    int ks;
    int l0;
    int conv;
    int mcpdg1;
    int mcpdg2;
    int mcppdg1;
    int mcppdg2;
    int mcpp1;
    int mcpp2;
  };
  VtxData _data;
};

#if 1
class MyAnalysis : public Algorithm {
  public:
    MyAnalysis(){}
    ~MyAnalysis() { _mcps = 0; }
    void init(Parameters* param);
    void process();
    void end();
  private:
    string _privtxname;
    string _secvtxname;
    string _v0vtxname;
    std::map<const MCParticle*,int> _mcpIndex;
    int _nEvt;

    //static const int Blist[]; 
    //static const int Clist[];
    //static const int Olist[];
    //static const int V0list[];
    static std::vector<int> Blist; 
    static std::vector<int> Clist;
    static std::vector<int> Olist;
    static std::vector<int> V0list;
    static const float aSmallNumber;

    class MCParticleExt
    {
      public:
        MCParticleExt() : _mcp(0), _ancestor(0), _descendant(0), _generation(0), _type(-1) 
                          {}
        MCParticleExt(const MCParticleExt& in) : _mcp(in._mcp), _ancestor(in._ancestor), _descendant(in._descendant),
                                        _generation(in._generation),
                                        _associatedTracks(in._associatedTracks), 
                                        _type(in._type) 
                                       {}
        ~MCParticleExt() { _mcp = 0; _ancestor = 0; _descendant = 0; } // Shouldn't deleted by this object.

        // getter
        const MCParticle*  getMCParticle() const {return _mcp;} 
        const MCParticle*  getAncestor()   const {return _ancestor;} 
        const MCParticle*  getDescendant() const {return _descendant;} 
        int          getGeneration() const {return _generation;} 
        int          getNAssociatedTracks() const {return _associatedTracks.size(); } 
        int          getType() const {return _type;} 
    //    TVector3     getRecPos() const {return _recPos;} 
        const Track* getAssociatedTrack(int itrk) const {return _associatedTracks[itrk]; } 

        // setter
        void setMCParticle(const MCParticle* in) { _mcp = in;} 
        void setAncestor(const MCParticle* in) { _ancestor = in;} 
        void setDescendant(const MCParticle* in) { _descendant = in;} 
        void setGeneration(int in) { _generation = in;} 
        void setType(int in) { _type = in;} 
//        void setRecPos(TVector3 in) { _recPos = in;} 

        void incrementGeneration() { _generation++;} 
        void decrementGeneration() { _generation--;} 
        void addAssociatedTrack(const Track* in) { _associatedTracks.push_back(in); } 

        //bool isAncestorOf(const MCParticle* in);
        //bool isDescendantOf(const MCParticle* in);

        static bool DescendingGeneration(const MCParticleExt& a, 
                                         const MCParticleExt& b) {
          return a.getGeneration() < b.getGeneration();
        }
        static bool AscendingGeneration(const MCParticleExt& a, 
                                        const MCParticleExt& b) {
          return b.getGeneration() < a.getGeneration();
        }
        static bool DescendingNAssociatedTracks(const MCParticleExt& a, 
                                                const MCParticleExt& b) {
          return a.getNAssociatedTracks() < b.getNAssociatedTracks();
        }
        static bool AscendingNAssociatedTracks(const MCParticleExt& a, 
                                               const MCParticleExt& b) {
//          int alpha = 0;
//          if (a.getNAssociatedTracks() == b.getNAssociatedTracks()) {
//            // if the numbers are same, then take the one that is nearer to MC vertex. 
//            std::cerr << "a.getNAssociatedTracks() = " << a.getNAssociatedTracks() << " b.getNAssociatedTracks() = " << b.getNAssociatedTracks() << std::endl;
//            std::cerr << "a.getAncestor() = " << a.getAncestor() << " b.getAncestor() = " << b.getAncestor() << std::endl;
//            std::cerr << "a.getRecPos().X() = " << a.getRecPos().X() << " b.getRecPos().X() = " << b.getRecPos().X() << std::endl;
            //std::cerr << "a.getAncestor()->getEndVertex().X() = " << a.getAncestor()->getEndVertex().X() << std::endl;
            //std::cerr << "(a.getRecPos() - a.getAncestor()->getEndVertex()).Mag() = " << (a.getRecPos() - a.getAncestor()->getEndVertex()).Mag() << std::endl;
            //std::cerr << "(b.getRecPos() - b.getAncestor()->getEndVertex()).Mag() = " << (b.getRecPos() - b.getAncestor()->getEndVertex()).Mag() << std::endl;
//            TVector3 mcPos_a;
//            if (a.getAncestor()->getPDG()==100000000 || a.getAncestor()->getPDG()==100000001) {
////std::cerr << "OK1" << std::endl;
//              mcPos_a = a.getAncestor()->getVertex();
//            } else mcPos_a = a.getAncestor()->getEndVertex();
//            TVector3 mcPos_b;
////std::cerr << "OK2" << std::endl;
//            if (b.getAncestor()->getPDG()==100000000 || b.getAncestor()->getPDG()==100000001) {
////std::cerr << "OK3" << std::endl;
//              mcPos_b = b.getAncestor()->getVertex();
//            } else mcPos_b = b.getAncestor()->getEndVertex();
////std::cerr << "OK4" << std::endl;
////std::cerr << "(a.getRecPos() - mcPos_a).Mag() = " << (a.getRecPos() - mcPos_a).Mag() << std::endl;
////std::cerr << "(b.getRecPos() - mcPos_b).Mag() = " << (b.getRecPos() - mcPos_b).Mag() << std::endl;
////bool test = false;
////if ((a.getRecPos() - mcPos_a).Mag() < (b.getRecPos() - mcPos_b).Mag()) test = true;
////std::cerr << "test = " << test << std::endl;
//            //return (a.getRecPos() - mcPos_a).Mag() < (b.getRecPos() - mcPos_b).Mag();
////            if ((b.getRecPos() - b.getAncestor()->getEndVertex()).Mag()<(a.getRecPos() - a.getAncestor()->getEndVertex()).Mag()) {
////              alpha++;
////            }
////          }  
//          //float alpha = -0.5;
//          //if ((b.getRecPos() - b.getAncestor()->getEndVertex()).Mag()<(a.getRecPos() - a.getAncestor()->getEndVertex()).Mag()) alpha = 0.5; 
          return b.getNAssociatedTracks() < a.getNAssociatedTracks();
//          float at = (a.getRecPos() - mcPos_a).Mag();
//          float bt = (b.getRecPos() - mcPos_b).Mag();
//std::cerr << "b.getNAssociatedTracks() + at/(at+bt) = " << b.getNAssociatedTracks() + at/(at+bt) << std::endl;
//std::cerr << "a.getNAssociatedTracks() + bt/(at+bt) = " << a.getNAssociatedTracks() + bt/(at+bt) << std::endl;
//          return b.getNAssociatedTracks() + at/(at+bt) < a.getNAssociatedTracks() + bt/(at+bt);
        }

       
      private:
        const MCParticle* _mcp;      
        const MCParticle* _ancestor;
        const MCParticle* _descendant;
        std::vector<const Track*>    _associatedTracks;
        int _generation;
        int _type; // origin type (0: beam bkg, 1: primary vertex, 2: B vertex, 3: C vertex, 4: O vertex)
//        TVector3 _recPos;
    };

    class MCVertex
    {
      public :
        MCVertex() : matchedvtxrec(0), ntrk(0), purity(0.) {}
        ~MCVertex() { 
           matchedvtxrec = 0; // Vertex object should NOT be deleted by this object.
        }
        MCParticleExt mcpe;
        const Vertex* matchedvtxrec;
        int ntrk;
        float purity;
    };

    // Data structure to fill tree.
    class EventData
    {
      public:
    	int nevt;
        float ipTruthx;
        float ipTruthy;
        float ipTruthz;
    };
    EventData _evdata;

    class VertexTrackData
    {
      public:
    	// mc informations
    	float  mcx_vtx;  // x of start vertex
    	float  mcy_vtx;  // y of start vertex
    	float  mcz_vtx;  // z of start vertex
    	float  mcx_evtx; // x of end vertex
    	float  mcy_evtx; // y of end vertex
    	float  mcz_evtx; // z of end vertex
    	double mce;
    	double mcp;
    	double mcpt;
    	double mcpx;
    	double mcpy;
    	double mcpz;
    	int    mcpdg;
    	double mccostheta; // Cos(theta)

    	// mc vtx (B-hadron/C-hadorn) informations
    	double mcve;  
    	double mcvp;
    	double mcvpt;
    	double mcvpx;
    	double mcvpy;
    	double mcvpz;
    	int    mcvpdg;
    	double mcvcostheta; // Cos(theta)

        int   isFromP;     // primary
        int   isFromB;     // B-hadrons
        int   isFromC;     // C-hadrons
        int   isFromO;     // other semi stables (e.g. gamma, tau,...)
        bool  isFromBmbkg; // beam background 
        float distancePvtxToSvtx;
    
    	// track informations
    	double rce;
    	double rcp;
    	double rcpt;
    	double rcpx;
    	double rcpy;
    	double rcpz;
    	double rccostheta; // Cos(theta)
	double d0;
	double d0err;
	double phi;
	double phierr;
	double omg;
	double omgerr;
	double z0;
	double z0err;
	double tanl;
	double tanlerr;

        bool isAssociatedToPvtx;
        bool isAssociatedToSvtx;
        bool isCorrectVertex;
        bool isCorrectDChain;

    };
    VertexTrackData _trkdata;

    class VertexData
    {
      public:
      float mcx;
      float mcy;
      float mcz;
      float mcpt;
      float mcpx;
      float mcpy;
      float mcpz;
      float mce;
      float mcp;
      float mccostheta;
      float mcmass;
      int   mcpdg;
      float ipTruthx;
      float ipTruthy;
      float ipTruthz;
      int   type;
      int   mcntrk;
      float weight;

      float rcx;
      float rcy;
      float rcz;
      float xerr;
      float yerr;
      float zerr;
      float rcpt;
      float rcpx;
      float rcpy;
      float rcpz;
      float rcp;
      float rce;
      float rccostheta;
      float rcmass;
      int   rcntrk;
      double chi2;

      double dl;   // residual in logitudinal to flight direction for secondary vertices.
      double dt1;  // residual in vertical to flight direction (dl X z) for secondary vertices.
      double dt2;  // residual in vertical to flight direction (dl X dt1) for secondary vertices.

      float distancePvtxToSvtx;
    };
    VertexData _vtxdata;

  public :
    bool              updateSemistableAncestorOf(MCParticleExt* p, std::vector<int>& semistables);
    bool              updateSemistableDescendantOf(MCParticleExt* p, std::vector<int>& semistables);
    const MCParticle* getParentIfSamePDG(const MCParticle* p);
    const MCParticle* getOriginMCParticle(const MCParticle* p);
    //TVector3          getIPTruth(const MCParticle* p);
    TVector3          getIPTruth();
    bool              isFromPrimaryVertex(const MCParticle* p);
    bool              isFromBeambackground(const MCParticle* p);
    void              checkVertices(std::vector<const Vertex*>& vtxs, 
                                    std::vector<MCVertex>& mcvtxs,
                                    std::map<const MCParticle*,MCParticleExt>& mcpMap,
                                    TTree* outtree);
    //void              setTracksOf(MCParticle* vtx, const MCParticleVec& mcps); 
    //void              setTracksOf(MCParticle* vtx); 
    MCParticle*       createVtxMCParticle(TVector3& pos); 
 private:
   const MCParticleVec* _mcps;
   TFile* _file;
   TTree* _trktree;
   TTree* _pvtxtree;
   TTree* _svtxtree;
   TTree* _v0vtxtree;
};
#endif
}

#endif

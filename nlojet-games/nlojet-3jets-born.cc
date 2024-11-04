/// This code computes the quantities we need from 2->2 Born-=level processes.
///


//------ DON'T TOUCH THIS PART! ------
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>
#include <bits/proc-hhc.h>  // for BSZ_extras
#include <nlo++-module_add.h>  // also for cmdline arguments 
#include "Numbers.hh"
// need this from NLOJET to reinitialise ran-numb seed.
#include "random.h"

#include "NLOHist.hh"
#include "helpers.hh"
#include <CmdLine.hh>

#include <cstdio>

//----- used namespaces -----
using namespace nlo;
using namespace std;

//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_hhc *, double&);
user_base_hhc * userfunc();

typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for hadron-hadron --> jets
	{"procindex", (void *) "hhc"},
  
	//   input function 
	{"inputfunc", (void *) inputfunc},
  
	//   phase space input function 
	{"psinput", (void *) psinput},
  
	//   user defined functions
	{"userfunc",  (void *) userfunc},
  
	//   module to generate the readable result
	{"main_module_add", (void *) module_add},
  
	//  end of the list
	{0, 0}
  };
}
//------ END OF THE DO-NOT-TOUCH-PART ------

//------ USER DEFINED PART STARTS HERE ------
#include <algorithm>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "pdf-genlha.h"
#pragma GCC diagnostic pop
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Boost.hh"

// try to see if we can get specialise our template for constructing
// an fj:PseudoJet from a HepMC::FourVector.
namespace fastjet {
  template<> PseudoJet::PseudoJet(const lorentzvector<double> & four_vector) {
    (*this) = PseudoJet(four_vector.X(),four_vector.Y(),four_vector.Z(),
			four_vector.T());
  }
}


using namespace fastjet;

const double sqrts=13600.0; // has to be declared as global

//------------------------------------------------------------------------
class UserHHC : public user1h_hhc{
public:
  //   init and user function
  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);
  
  // should overload the version in basic_user
  void operations_at_the_end_of_event(); 

private:
  //   pdf
  //pdf_cteq6 pdf;
  pdf_and_coupling_hhc * pdfp;
  
  //  conversion form weight_hhc to double
  weight_conversion<weight_hhc, double> conversion;

  vector<double> murs;
  vector<double> mufs;

  // dijet selection
  const double dijet_R      = 0.4;
  const double dijet_ymax   = 1.7;
  const double dijet_ptmin1 = 700.0;
  const double dijet_ptmin2 = 150.0;
  const double dijet_ptmax2 = 200.0;
  const double dijet_dRmin  = 1.0;
  const double dijet_dRmax  = 1.2;

  // grooming selection
  const double groom_R     = 1.2;
  const double groom_ymax  = 1.7;
  const double groom_ptmin = 1000.0;
  const double groom_zmin  = 0.1;
  const double groom_zmax  = 0.2;
  const double groom_dRmin = 0.8;


  // binning info
  vector<NLOHist> histograms; // each histogram corresponds to a slice in rapidity
  enum HistIndices{
    dijet_tot = 0,
    dijet_g   = 1,
    groom_tot = 2,
    groom_g   = 3,
  };
  const unsigned int nhist = 4;
 
  unsigned int iev;
  int __save_after;
  string __ohistfile;
  ostringstream _header;

  unsigned int _pdf_channel(int p1, int p2) const;

};


//----- defines the module to sum up the results of the different runs -----
module_add_type module_add = main_module_add<basic_user_result<user1h_hhc::distbook_type> >;

user_base_hhc * userfunc() {
  return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd){
  // number of jets
  nj = 3U;

  // number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_hhc *ps, double& s){
  // total c.m. energy square
  s = sqrts*sqrts;
  
  // You can use your own phase generator. 
  // Here we use the default.
  ps = 0;
} 


void UserHHC::initfunc(unsigned int){
  // command-line parsing
  //---------------------
  int argc;
  char **argv;
  findargs(&argc, &argv);
  CmdLine cmd(argc, argv);
  
  // jet radius and min pt cut
  //

  // randopm seed
  unsigned int iseq = (unsigned int)(cmd.value("-iseq",1));
  __save_after = int(cmd.value("--save-after", 10000));
  string pdfname = cmd.value<string>("--pdf", "PDF4LHC21_mc");
  int    pdf_mem = int(cmd.value("--pdfmem", 0));


  // figure out what level we're dealing with
  string level = cmd.value<string>("-c");
  
    // create the header
  _header << "# Ran: " << cmd.command_line() << endl;
  _header << "#" << endl;
  _header << "# Parameters:" << endl;
  _header << "#   sqrts      = " << sqrts   << endl;
  _header << "#" << endl;
  _header << "#   level      = " << level << endl;
  _header << "#   iseq       = " << iseq << endl;
  _header << "#   save_after = " << __save_after << endl;
  _header << "#----------------------------------------" << endl;

  // deal with pdfs...
  //string pdfname = "PDF4LHC21_mc";
  string fullname = pdfname;
  //int    pdf_mem  = 0;
  pdfp = new pdf_genlha_pp(fullname,pdf_mem);

  // output/parts filename
  ostringstream oss;
  if (pdf_mem>0){
    oss << "output/parts/3j-" << pdfname << "_" << pdf_mem << "-" << level << "-iseq" << iseq << ".ohist";
  } else {
    oss << "output/parts/3j-" << pdfname << "-" << level << "-iseq" << iseq << ".ohist";
  }
  __ohistfile = oss.str();

  cout << "results will be saved to: " << endl << "  " << __ohistfile << endl;


  // decide the scales we'll use
  murs.push_back(1.0);  mufs.push_back(1.0);
  murs.push_back(0.5);  mufs.push_back(0.5);
  murs.push_back(0.5);  mufs.push_back(1.0);
  murs.push_back(1.0);  mufs.push_back(0.5);
  murs.push_back(1.0);  mufs.push_back(2.0);
  murs.push_back(2.0);  mufs.push_back(1.0);
  murs.push_back(2.0);  mufs.push_back(2.0);
 
  // binning details
  //---------------- 
  histograms.resize(nhist*mufs.size()*3);
  for (unsigned int iscale=0; iscale<mufs.size(); ++iscale){
    ostringstream oss;
    oss << "-mur_" << murs[iscale] << "-muf_" << mufs[iscale];
    
    histograms[(nhist*iscale+dijet_tot)*3+0].declare("sigma_dijet_v_ptmin"     +oss.str(), -200, -100, 10);
    histograms[(nhist*iscale+dijet_tot)*3+1].declare("sigma_dijet_v_ptmax"     +oss.str(),  150,  300, 15);
    histograms[(nhist*iscale+dijet_tot)*3+2].declare("sigma_dijet_v_drmin"     +oss.str(), -1.2, -0.0, 12);
    histograms[(nhist*iscale+dijet_g  )*3+0].declare("sigma_dijet_g_v_ptmin"   +oss.str(), -200, -100, 10);
    histograms[(nhist*iscale+dijet_g  )*3+1].declare("sigma_dijet_g_v_ptmax"   +oss.str(),  150,  300, 15);
    histograms[(nhist*iscale+dijet_g  )*3+2].declare("sigma_dijet_g_v_drmin"   +oss.str(), -1.2, -0.0, 12);

    histograms[(nhist*iscale+groom_tot)*3+0].declare("sigma_groom_v_zmin"      +oss.str(), -0.2, -0.05, 15);
    histograms[(nhist*iscale+groom_tot)*3+1].declare("sigma_groom_v_zmax"      +oss.str(),  0.1,  0.5 , 20);
    histograms[(nhist*iscale+groom_tot)*3+2].declare("sigma_groom_v_drmin"     +oss.str(), -1.2, -0.0 , 12);
    histograms[(nhist*iscale+groom_g  )*3+0].declare("sigma_groom_g_v_zmin"    +oss.str(), -0.2, -0.05, 15);
    histograms[(nhist*iscale+groom_g  )*3+1].declare("sigma_groom_g_v_zmax"    +oss.str(),  0.1,  0.5 , 20);
    histograms[(nhist*iscale+groom_g  )*3+2].declare("sigma_groom_g_v_drmin"   +oss.str(), -1.2, -0.0 , 12);

  }

  cout << "End-user initilisation done" << endl; 

}

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp){
  // first make a quick test that ht is large enough
  //
  // We need at least one jet with pt>600
  // so the scalar sum should be at least 
  // and the 2nd pt has to be at least 2*600
  double htsum = 0.0;
  for (int i = 1; i <= p.upper(); i++) {
    htsum += p[i].perp();
  }
  if (htsum<2*dijet_ptmin1) return;

  assert(p.upper() == 3); //< hardcoded for 3 particles in the final state

  //--------------------------------------------------
  // transfer event into our own format
  // nlojet partons 1, 2, 3
  // FJ     partons 0, 1, 2
  // fj indices     1, 2, 3
  vector<PseudoJet> partons(p.upper());
  for (int i = 1; i <= p.upper(); i++) {
    partons[i-1] = p[i];
    partons[i-1].set_user_index(i);

    // give up on events where momenta have NaN
    if ( isbad(partons[i-1].E()) ) return;
  }
  // eliminate (z-momentum balancing) final parton with zero pt (redundant).
  if (partons.back().pt2() == 0.0) partons.pop_back();

  //----------------------------------------------------------------------
  // helpers to compute the weights
  //----------------------------------------------------------------------
  // scales
  double ref_scaleR2 = htsum * htsum; // is this what we want?
  double ref_scaleF2 = htsum * htsum;

  // NLOJet++ weights
  vector<weight_hhc> nlo_weight_vectors;
  vector<double> nlo_weights;

  // weight fractions for each flavour channel  
  // 
  // In a nutshell, we have 14 flavour channel mapped onto 7 PDF
  // channels.  We need to know the relative weight of each of the
  // flavour channel within its associated PDF channel
  //
  // This info is stored in flv_weights
  vector<double> pdf_weights(7, 0.0);
  vector<double> flv_weights(14, 0.0);
  bool weights_computed = false;

  auto compute_flavoured_weights = [&](){
    // NLOJet weights
    for (unsigned iscale = 0; iscale < murs.size(); iscale++){
      double rscale2 = ref_scaleR2 * murs[iscale] * murs[iscale];
      double fscale2 = ref_scaleF2 * mufs[iscale] * mufs[iscale];
      nlo_weight_vectors.push_back(amp(*pdfp, rscale2, fscale2, 389385.730));
      nlo_weights.push_back(conversion(nlo_weight_vectors.back()));
    }

    // flavour info
    for (unsigned ichan = 0; ichan < BSZ_flav.ndecomp(); ichan++){
      BSZ_event_flav qrs_flav = BSZ_flav.flav(ichan);
      unsigned int pdf_channel = _pdf_channel(qrs_flav[-1],qrs_flav[0]);
      pdf_weights[pdf_channel] += BSZ_flav.weight(ichan);
      flv_weights[ichan]       += BSZ_flav.weight(ichan);
    }
    for (unsigned ichan = 0; ichan < BSZ_flav.ndecomp(); ichan++){
      BSZ_event_flav qrs_flav = BSZ_flav.flav(ichan);
      unsigned int pdf_channel = _pdf_channel(qrs_flav[-1],qrs_flav[0]);
      flv_weights[ichan] /= pdf_weights[pdf_channel];
    }
    weights_computed = true;
  };

  auto compute_gluon_weights = [&](unsigned int nlojet_index) -> vector<double>{
    if (!weights_computed) compute_flavoured_weights();
    
    // then we loop over the flavour channels
    vector<weight_hhc> nlo_weight_gluon_vectors(murs.size());
    vector<double> nlo_weights_gluon(murs.size(), 0.0);
    for (unsigned ichan = 0; ichan < BSZ_flav.ndecomp(); ichan++){
      // flavour info for this channel
      BSZ_event_flav qrs_flav = BSZ_flav.flav(ichan);
      
      // only gluons are needed
      if (qrs_flav[nlojet_index]!=0) continue;
      
      // pdf channel from incoming flavours
      unsigned int pdf_channel = _pdf_channel(qrs_flav[-1],qrs_flav[0]);
      for (unsigned iscale = 0; iscale < murs.size(); iscale++){
        nlo_weight_gluon_vectors[iscale][pdf_channel] += nlo_weight_vectors[iscale][pdf_channel]*flv_weights[ichan];
      }
    }
    
    // convolute w PDFs
    for (unsigned iscale = 0; iscale < murs.size(); iscale++){
      nlo_weights_gluon[iscale] = conversion(nlo_weight_gluon_vectors[iscale]);
    }
    return nlo_weights_gluon;
  };
  
  
  //----------------------------------------------------------------------
  // dijet selection
  //----------------------------------------------------------------------
  // cluster the event
  JetDefinition dijet_jet_def(antikt_algorithm, dijet_R);
  Selector sel_dijets = SelectorAbsRapMax(dijet_ymax);
  vector<PseudoJet> dijet_jets = sel_dijets(dijet_jet_def(partons));

  // browse all pairs
  for (unsigned int i1=0; i1<dijet_jets.size(); ++i1){
    for (unsigned int i2=i1+1; i2<dijet_jets.size(); ++i2){
      // impose common cuts
      //..............................
      // pt cut on the leading jet
      if (dijet_jets[i1].pt2()<dijet_ptmin1*dijet_ptmin1) continue;
      // delta R between the jets
      double dR2 = dijet_jets[i1].squared_distance(dijet_jets[i2]);
      if (dR2>dijet_dRmax*dijet_dRmax) continue;

      bool pass_ptmin2 = (dijet_jets[i2].pt2()>dijet_ptmin2*dijet_ptmin2);
      bool pass_ptmax2 = (dijet_jets[i2].pt2()<dijet_ptmax2*dijet_ptmax2);
      bool pass_drmin  = (dR2>dijet_dRmin*dijet_dRmin);

      unsigned int nlojet_index = dijet_jets[i2].user_index();
      
      // variation with ptmin2
      //..............................
      if (pass_ptmax2 && pass_drmin){
        vector<double> nlo_weights_gluon = compute_gluon_weights(nlojet_index);
        double pt2 = dijet_jets[i2].pt();
        for (unsigned iscale = 0; iscale < murs.size(); iscale++){
          histograms[(nhist*iscale+dijet_tot)*3+0].add_entry(-pt2, nlo_weights      [iscale]);
          histograms[(nhist*iscale+dijet_g  )*3+0].add_entry(-pt2, nlo_weights_gluon[iscale]);
        }        
      }
      
      // variation with ptmax2
      //..............................
      if (pass_ptmin2 && pass_drmin){
        vector<double> nlo_weights_gluon = compute_gluon_weights(nlojet_index);
        double pt2 = dijet_jets[i2].pt();
        for (unsigned iscale = 0; iscale < murs.size(); iscale++){
          histograms[(nhist*iscale+dijet_tot)*3+1].add_entry(pt2, nlo_weights      [iscale]);
          histograms[(nhist*iscale+dijet_g  )*3+1].add_entry(pt2, nlo_weights_gluon[iscale]);
        }
      }

      // variation with drmin
      //..............................
      if (pass_ptmin2 && pass_ptmax2){
        vector<double> nlo_weights_gluon = compute_gluon_weights(nlojet_index);
        double dR = sqrt(dR2);
        for (unsigned iscale = 0; iscale < murs.size(); iscale++){
          histograms[(nhist*iscale+dijet_tot)*3+2].add_entry(-dR, nlo_weights      [iscale]);
          histograms[(nhist*iscale+dijet_g  )*3+2].add_entry(-dR, nlo_weights_gluon[iscale]);
        }
      }

      
    }
  }      
  
  //----------------------------------------------------------------------
  // groomed selection
  //----------------------------------------------------------------------
  // cluster the event
  JetDefinition groom_jet_def(antikt_algorithm, groom_R);
  Selector sel_groom = SelectorAbsRapMax(groom_ymax) * SelectorPtMin(groom_ptmin);
  vector<PseudoJet> groom_jets = sel_dijets(groom_jet_def(partons));

  // browse all jets
  for (const auto & jet : groom_jets){
    // impose common cuts
    //..............................
      
    // we need at least 2 constituents
    const auto constits = sorted_by_pt(jet.constituents());
    if (constits.size()<2) continue;

    // zcut
    double pt1 = constits[0].pt();
    double pt2 = constits[1].pt();
    double z = pt2/(pt1+pt2);

    bool pass_zmin = (z>groom_zmin);
    bool pass_zmax = (z<groom_zmax);
    
    // rg cut
    double dR2 = constits[0].squared_distance(constits[1]);
    bool pass_drmin = (dR2>groom_dRmin*groom_dRmin);

    unsigned int nlojet_index = constits[1].user_index();

    // variation with ptmin2
    //..............................
    if (pass_zmax && pass_drmin){
      vector<double> nlo_weights_gluon = compute_gluon_weights(nlojet_index);
      for (unsigned iscale = 0; iscale < murs.size(); iscale++){
        histograms[(nhist*iscale+groom_tot)*3+0].add_entry(-z, nlo_weights      [iscale]);
        histograms[(nhist*iscale+groom_g  )*3+0].add_entry(-z, nlo_weights_gluon[iscale]);
      }        
    }
    
    // variation with ptmax2
    //..............................
    if (pass_zmin && pass_drmin){
      vector<double> nlo_weights_gluon = compute_gluon_weights(nlojet_index);
      for (unsigned iscale = 0; iscale < murs.size(); iscale++){
        histograms[(nhist*iscale+groom_tot)*3+1].add_entry(z, nlo_weights      [iscale]);
        histograms[(nhist*iscale+groom_g  )*3+1].add_entry(z, nlo_weights_gluon[iscale]);
      }
    }
    
    // variation with drmin
    //..............................
    if (pass_zmin && pass_zmax){
      vector<double> nlo_weights_gluon = compute_gluon_weights(nlojet_index);
      double dR = sqrt(dR2);
      for (unsigned iscale = 0; iscale < murs.size(); iscale++){
        histograms[(nhist*iscale+groom_tot)*3+2].add_entry(-dR, nlo_weights      [iscale]);
        histograms[(nhist*iscale+groom_g  )*3+2].add_entry(-dR, nlo_weights_gluon[iscale]);
      }
    }
  }    

}


void UserHHC::operations_at_the_end_of_event(){
  
  // deal with output/parts
  if (histograms.size() > 0) {

    // collate histograms
    for (auto &hist : histograms)  hist.collate_event();
    
    // figure out number of events
    double nev = histograms.begin()->nev();
    double nev_bad = histograms.begin()->nev_bad();
    if (abs(nev - __save_after * int(nev/__save_after+0.5)) < 0.1) {
      
      ofstream fstr( __ohistfile.c_str() );

      fstr << _header.str();

      fstr << "# nev = "     << nev     << endl;
      fstr << "# nev_bad = " << nev_bad << endl;
      for (auto &hist : histograms){
	hist.write(fstr, 0);
	fstr << endl << endl;
      }

    }
  }

  // remeber to run NLOjet's official end of event operations
  user1h_hhc::operations_at_the_end_of_event();
}

unsigned int UserHHC::_pdf_channel(int p1, int p2) const{
  // 0: gg
  // 1: qg
  // 2: gq
  // 3: qr
  // 4: qq
  // 5: qqb
  // 6: qrb
  if (p1==0){
    if      (p2==0){ return 0; } // gg
    else if (p2==1){ return 2; } // gq
    else { 
      cerr << "p1=0, p2=" << p2  << " is not a proper pdf flavour channel" << endl;
      exit(1);
    }
  } else if (p1==1){
    switch (p2){
    case  0: return 1; break; // qg
    case  1: return 4; break; // qq
    case -1: return 5; break; // qqb
    case  2: return 3; break; // qr
    case -2: return 6; break; // qrb
    default: 
      cerr << "p1=1, p2=" << p2  << " is not a proper pdf flavour channel" << endl;
      exit(1);
    };
  } else {
    cerr << "p1=" << p1 << " is not a proper pdf flavour channel" << endl;
    exit(1);
  }
}


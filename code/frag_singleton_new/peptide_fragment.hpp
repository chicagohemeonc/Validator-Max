#ifndef __peptide_fragment_hpp_
#define __peptide_fragment_hpp_

#include <string>
#include <list>
using namespace std;

extern "C" double compute_pI(char *seq, unsigned long seq_length, int charge_increment);

#define SAMPLEPEPTIDE "DIGSESTEDQAMEDIK"
#define MAX_SEQUENCE 100000
#define MAX_CHARGE   10
#define PROTON_MASS 1.00727646688

struct ProteaseStruct
{
   char szName[96];
   char szCut[24];
   char szNoCut[24];
   int  iSense;    /* 0=n-term, 1=c-term */
} *pProtease, *pTmp;


typedef std::list<double> ListD;


class PeptideFragment {

public:
  PeptideFragment ( int massType );
  void analyze( char * peptide, int charge, char *modification );
  string sequence ( void ) const;
  virtual ListD a_ions ( void );
  virtual ListD b_ions ( void );
  virtual ListD y_ions ( void );
  virtual ListD peptide_mass ( void );
  virtual string composition ( void );
  virtual double pI ( void );
  ListD a;
  ListD b;
  ListD c;
  ListD x;
  ListD y;
  ListD z;
  ListD zdot;
  ListD mass;

private:
  /* parameter values */
  char _szInputSequence[MAX_SEQUENCE];
  unsigned int _iLenPeptide;
  unsigned _iMassType; /*
                   * bMassType=0 for average masses
                   * bMassType=1 for monoisotopic masses
                   */
  int _iCharge;          /* ion charge */
  string _modification;  /* specifies which residues are modified and how */

  /* work variables */
  double pdMassAA[128];
  double pdPositionMod[MAX_SEQUENCE];
  double dPepMass;

  int piCompC[128];
  int piCompH[128];
  int piCompN[128];
  int piCompO[128];
  int piCompS[128];

  /* private methods */
  void _initialize_mass ( void );
  void _init_comp ( void );
};

#endif
